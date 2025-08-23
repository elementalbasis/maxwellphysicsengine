using DifferentialEquations
using LinearAlgebra
using Printf
using StaticArrays
using UUIDs
import Base: @kwdef


###############################################################################
# SYSTEM BASICS
###############################################################################

O = [0.0, 0.0, 0.0]
X = [1.0, 0.0, 0.0]
Y = [0.0, 1.0, 0.0]
Z = [0.0, 0.0, 1.0]

abstract type Entity end
abstract type Body <: Entity end
abstract type Force <: Entity end

@kwdef mutable struct System
	forces::Vector{Force} = []
	bodies::Vector{Body} = []
	entities::Vector{Entity} = []
	state::Vector{Float64} = []
	time::Float64 = 0.0
	index_map::Dict{Entity, UnitRange{Int}} = Dict()
end

function get_state(system::System, entity::Entity,
		system_state::Vector{Float64})
	index_range = system.index_map[entity]
	return system_state[index_range]
end
get_state(system::System,
	  entity::Entity) = get_state(system, entity, system.state)

function set_state!(system::System, entity::Entity,
		entity_state::Vector{Float64})
	index_range = system.index_map[entity]
	system.state[index_range] = entity_state
end

function register!(system::System, entity::Entity)
	push!(system.entities, entity)
	entity isa Force && push!(system.forces, entity)
	entity isa Body && push!(system.bodies, entity)

	n = entity_state_size(entity)
	n == 0 && return

	index_start = length(system.state) + 1
	index_stop = index_start + n - 1
	index_range = index_start : index_stop

	append!(system.state, zeros(n))
	system.index_map[entity] = index_range
end



###############################################################################
# Particles
###############################################################################

@kwdef struct Particle <: Body
	uuid::UUID = uuid4()
	mass::Float64 = 1.0
end
entity_state_size(particle::Particle) = 6

function create_particle!(system::System;
		position::Vector{Float64} = O,
		velocity::Vector{Float64} = O,
		mass::Float64 = 1.0,
		)
	particle = Particle(mass = mass)
	register!(system, particle)

	particle_state = [position; velocity]
	set_state!(system, particle, particle_state)

	return particle
end

function get_position(system::System, particle::Particle,
		system_state::Vector{Float64})
	particle_state = get_state(system, particle, system_state)
	return particle_state[1:3]
end
get_position(system::System,
	     particle::Particle) = get_position(system, particle, system.state)

function get_velocity(system::System, particle::Particle,
		system_state::Vector{Float64})
	particle_state = get_state(system, particle, system_state)
	return particle_state[4:6]
end
get_velocity(system::System,
	     particle::Particle) = get_velocity(system, particle, system.state)

function set_position!(system::System, particle::Particle, position::Vector{Float64})
	particle_state = get_state(system, particle)
	particle_state[1:3] = position
	set_state!(system, particle, particle_state)
end

function set_velocity!(system::System, particle::Particle, velocity::Vector{Float64})
	particle_state = get_state(system, particle)
	particle_state[4:6] = velocity
	set_state!(system, particle, particle_state)
end



###############################################################################
# Forces
###############################################################################

# General

compute_force(system::System, force::Force,
	      body::Body) = compute_force(system, force, body, system.state)



# Uniform Gravity

@kwdef struct UniformGravity <: Force
	g = 9.8
	direction = -Z
	targets = :all
end

function compute_force(system::System, force::UniformGravity, body::Body,
		system_state::Vector{Float64})
	return body.mass * force.g * force.direction
end



# Linear Drag

@kwdef struct LinearDrag <: Force
	b = 1
	targets = :all
end

function compute_force(system::System, force::LinearDrag, body::Body,
		system_state::Vector{Float64})
	v = get_velocity(system, body, system_state)
	return -v * force.b
end



# Spring

@kwdef struct Spring <: Force
	k = 1
	length = 0
	targets::Tuple{Particle, Particle}
end

function compute_force(system::System, force::Spring, particle::Particle,
		system_state::Vector{Float64})
	particle in force.targets || return O

	pa = force.targets[1]
	pb = force.targets[2]
	ra = get_position(system, pa, system_state)
	rb = get_position(system, pb, system_state)
	u = rb - ra
	x = norm(u)
	n = (x == 0) ? O : normalize(u)
	f = force.k * (x - force.length) * n

	return Dict(pa => f, pb => -f)[particle]
end



# Modulated Spring

@kwdef struct ModulatedSpring <: Force
	k_max = 1
	sensitivity = 1
	length = 1
	targets::Tuple{Particle, Particle}
end
entity_state_size(force::ModulatedSpring) = 1

function get_modulated_spring_activation(system::System, force::ModulatedSpring,
		system_state::Vector{Float64})
	entity_state = get_state(system, force, system_state)
	return entity_state[1]
end
get_modulated_spring_activation(
	system::System, force::ModulatedSpring) = get_state(
		system, force, system.state)

function set_modulated_spring_activation!(system::System,
		force::ModulatedSpring, activation::Float64)
	set_state!(system, force, [activation])
end

function compute_force(system::System, force::ModulatedSpring,
		particle::Particle, system_state::Vector{Float64})
	particle in force.targets || return O

	pa = force.targets[1]
	pb = force.targets[2]
	ra = get_position(system, pa, system_state)
	rb = get_position(system, pb, system_state)

	u = rb - ra
	x = norm(u)
	n = (x == 0) ? O : normalize(u)

	k = get_modulated_spring_activation(system, force, system_state)
	f = k * x

	return Dict(pa => f, pb => -f)[particle]
end

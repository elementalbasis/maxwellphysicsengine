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
	index_map::Dict{UUID, UnitRange{Int}} = Dict()
end

function get_state(system::System, entity::Entity,
		system_state::Vector{Float64})
	index_range = system.index_map[entity.uuid]
	return system_state[index_range]
end
get_state(system::System,
	  entity::Entity) = get_state(system, entity, system.state)

function set_state!(system::System, entity::Entity,
		entity_state::Vector{Float64})
	index_range = system.index_map[entity.uuid]
	system.state[index_range] = entity_state
end

function register_entity!(system::System, entity::Entity)
	push!(system.entities, entity)
	if entity isa Force
		push!(system.forces, entity)
	elseif entity isa Body
		push!(system.bodies, entity)
	end

	n = entity_state_size(entity)
	if n == 0
		return
	end

	index_start = length(system.state) + 1
	index_stop = index_start + n - 1
	index_range = index_start : index_stop

	append!(system.state, zeros(n))
	system.index_map[entity.uuid] = index_range
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
	register_entity!(system, particle)

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
# Force
###############################################################################

# General

compute_force(system::System,
	      force::Force,
	      particle::Particle) = compute_force(system, force, particle,
						  system.state)



# Uniform Gravity

@kwdef struct UniformGravity <: Force
	g = 9.8
	direction = -Z
end

acts_on(system::System, force::UniformGravity, particle::Particle) = true
function compute_force(system::System, force::UniformGravity, particle::Particle,
		system_state::Vector{Float64})
	return particle.mass * force.g * force.direction
end



# Linear Drag

@kwdef struct LinearDrag <: Force
	b = 1
end

acts_on(system::System, force::LinearDrag, particle::Particle) = true
function compute_force(system::System, force::LinearDrag, particle::Particle,
		system_state::Vector{Float64})
	v = get_velocity(system, particle, system_state)
	return -v * force.b
end

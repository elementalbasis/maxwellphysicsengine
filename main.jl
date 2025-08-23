using DifferentialEquations
using LinearAlgebra
using Printf
using StaticArrays
using StatsFuns
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

function get_state(system::System, entity::Entity;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	if system_state == nothing
		system_state = system.state
	end
	index_range = system.index_map[entity]
	return system_state[index_range]
end

function set_state!(system::System, entity::Entity,
		entity_state::Vector{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	if system_state == nothing
		system_state = system.state
	end
	index_range = system.index_map[entity]
	system_state[index_range] = entity_state
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

function register!(system::System, entities::Entity...)
	for e in entities
		register!(system, e)
	end
end

entity_state_size(entity::Entity) = 0



###############################################################################
# Particles
###############################################################################

@kwdef struct Particle <: Body
	uuid::UUID = uuid4()
	mass::Float64 = 1.0
	is_stationary = false
end
entity_state_size(particle::Particle) = 6

function create_particle!(system::System;
		position::Vector{Float64} = O,
		velocity::Vector{Float64} = O,
		kwargs...
		)
	particle = Particle(;kwargs...)
	register!(system, particle)

	particle_state = [position; velocity]
	set_state!(system, particle, particle_state)

	return particle
end

function get_position(system::System, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	particle_state = get_state(system, particle, system_state = system_state)
	return particle_state[1:3]
end

function get_velocity(system::System, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	particle_state = get_state(system, particle, system_state = system_state)
	return particle_state[4:6]
end

function set_position!(system::System, particle::Particle, position::Vector{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	particle_state = get_state(system, particle, system_state = system_state)
	particle_state[1:3] = position
	set_state!(system, particle, particle_state, system_state = system_state)
end

function set_velocity!(system::System, particle::Particle, velocity::Vector{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	particle_state = get_state(system, particle, system_state = system_state)
	particle_state[4:6] = velocity
	set_state!(system, particle, particle_state, system_state = system_state)
end

function get_acceleration(system::System, body::Body;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	body.is_stationary && return O
	F = O
	for force in system.forces
		F += compute_force(system, force, body, system_state = system_state)
	end

	a = F / body.mass
	return a
end

function get_state_flow(system::System, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	v = get_velocity(system, particle, system_state = system_state)
	a = get_acceleration(system, particle, system_state = system_state)
	return [v; a]
end

###############################################################################
# Forces
###############################################################################



# Uniform Gravity

@kwdef struct UniformGravity <: Force
	g = 9.8
	direction = -Z
	targets = :all
end

function compute_force(system::System, force::UniformGravity, body::Body;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	return body.mass * force.g * force.direction
end



# Linear Drag

@kwdef struct LinearDrag <: Force
	b = 1
	targets = :all
end

function compute_force(system::System, force::LinearDrag, body::Body;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	v = get_velocity(system, body, system_state = system_state)
	return -v * force.b
end



# Spring

@kwdef struct Spring <: Force
	k = 1
	length = 0
	targets::Tuple{Particle, Particle}
end

function compute_force(system::System, force::Spring, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	particle in force.targets || return O

	pa = force.targets[1]
	pb = force.targets[2]
	ra = get_position(system, pa, system_state = system_state)
	rb = get_position(system, pb, system_state = system_state)
	u = rb - ra
	x = norm(u)
	n = (x == 0) ? O : normalize(u)
	f = force.k * (x - force.length) * n

	return Dict(pa => f, pb => -f)[particle]
end



# Modulated Spring

@kwdef struct ModulatedSpring <: Force
	k_max = 10
	sensitivity = 0.01
	length = 1
	targets::Tuple{Particle, Particle}
end
entity_state_size(force::ModulatedSpring) = 1

function get_modulated_spring_activation(system::System, force::ModulatedSpring;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	entity_state = get_state(system, force, system_state = system_state)
	return entity_state[1]
end

function set_modulated_spring_activation!(system::System,
		force::ModulatedSpring, activation::Float64;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	set_state!(system, force, [activation], system_state = system_state)
end

function get_state_flow(system::System, force::ModulatedSpring;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	pa = force.targets[1]
	pb = force.targets[2]
	ra = get_position(system, pa, system_state = system_state)
	rb = get_position(system, pb, system_state = system_state)

	u = rb - ra
	x = norm(u)

	return [force.sensitivity * (x - force.length)]
end


function compute_force(system::System, force::ModulatedSpring,
		particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	particle in force.targets || return O

	pa = force.targets[1]
	pb = force.targets[2]
	ra = get_position(system, pa, system_state = system_state)
	rb = get_position(system, pb, system_state = system_state)

	u = rb - ra
	x = norm(u)
	n = (x == 0) ? O : normalize(u)

	q = get_modulated_spring_activation(system, force,
					    system_state = system_state)
	k = force.k_max * logistic(q)

	f = k * x * n

	return Dict(pa => f, pb => -f)[particle]
end




###############################################################################
# Solver
###############################################################################

function get_state_flow(system::System;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	n = length(system.state)
	system_state_flow = zeros(n)
	for entity in keys(system.index_map)
		entity_state_flow = get_state_flow(system, entity,
						   system_state = system_state)
		set_state!(system, entity, entity_state_flow,
			   system_state = system_state_flow)
	end
	return system_state_flow
end

function ODE_function(system_state::Vector{Float64}, system::System, t::Float64)
	return get_state_flow(system, system_state = system_state)
end

function update!(system::System, time_increment)
	tspan = (system.time, system.time + time_increment)
	prob = ODEProblem(ODE_function, system.state, tspan, system)
	sol = solve(prob)
	system.time += time_increment
	system.state = sol(system.time)

	return sol
end

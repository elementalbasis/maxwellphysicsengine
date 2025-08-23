#using DifferentialEquations
#using LinearAlgebra
#using Printf
#using StaticArrays
#using StatsFuns
using UUIDs
import Base: @kwdef



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

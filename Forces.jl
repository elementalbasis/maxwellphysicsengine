#using DifferentialEquations
using LinearAlgebra
#using Printf
#using StaticArrays
using StatsFuns
#using UUIDs
import Base: @kwdef



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

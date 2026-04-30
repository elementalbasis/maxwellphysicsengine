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



# Particle Contact

@kwdef struct ParticleContact <: Force
	k = 1
	chunks::Chunks
end


function compute_force(system::System, force::ParticleContact, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	ensure_updated!(force.chunks, system, system_state = system_state)

	ra = get_position(system, particle, system_state = system_state)
	Ra = particle.radius
	F = O
	neighborhood = nearby_particles(force.chunks, particle, system, system_state = system_state)
	for pb in neighborhood
		rb = get_position(system, pb, system_state = system_state)
		Rb = particle.radius
		R = Ra + Rb
		u = rb - ra
		d = norm(u)
		n = (d == 0) ? O : normalize(u)
		if d < R
			F += force.k * (-n) * (R - d)
		end
	end
	return F
end



# Wall Contact

# This is a plane that has an allowed side and a prohibited side
@kwdef struct WallContact <: Force
	k = 1
	n = O # A unit vector that points in the direction of the allowed side
	p = O # A point on the plane
end

function compute_force(system::System, force::WallContact, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	ra = get_position(system, particle, system_state = system_state)
	Ra = particle.radius
	s = dot(force.n, (ra - force.p))
	if (s > 0)
		return O
	else
		return force.n * (-s) * force.k
	end
end



# Modulated Wall Contact

# This is a plane that has an allowed side and a prohibited side
@kwdef struct ModulatedWallContact <: Force
	k = 1
	b_max = 0.001
	sensitivity = 0.0001
	energy_target = 0.01
	n = O # A unit vector that points in the direction of the allowed side
	p = O # A point on the plane
end
entity_state_size(force::ModulatedWallContact) = 1

function get_modulated_wall_contact_activation(system::System, force::ModulatedWallContact;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	entity_state = get_state(system, force, system_state = system_state)
	return entity_state[1]
end

function set_modulated_wall_contact_activation!(system::System,
		force::ModulatedWallContact, activation::Float64;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	set_state!(system, force, [activation], system_state = system_state)
end

function get_state_flow(system::System, force::ModulatedWallContact;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	K = get_average_kinetic_energy(system, system_state = system_state)
	return [force.sensitivity * (K - force.energy_target)]
end

function compute_force(system::System, force::ModulatedWallContact, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	ra = get_position(system, particle, system_state = system_state)
	va = get_velocity(system, particle, system_state = system_state)
	Ra = particle.radius
	s = dot(force.n, (ra - force.p))
	q = get_modulated_wall_contact_activation(system, force, system_state = system_state)
	b = force.b_max * logistic(q)
	if (s > 0)
		return O
	else
		return force.n * (-s) * force.k - b * va
	end
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

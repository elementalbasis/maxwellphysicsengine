#using DifferentialEquations
using LinearAlgebra
#using Printf
#using StaticArrays
#using StatsFuns
using Rotations
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



###############################################################################
# Rigid Bodies
###############################################################################


@kwdef struct RigidBody <: Body
	uuid::UUID = uuid4()
	mass::Float64 = 1.0
	Jb::Matrix{Float64} = Matrix{Float64}(I, 3, 3)
	invJb::Matrix{Float64} = Matrix{Float64}(I, 3, 3)
	is_stationary = false
end
entity_state_size(rb::RigidBody) = 13

function create_rigid_body!(system::System;
		position::Vector{Float64} = O,
		velocity::Vector{Float64} = O,
		orientation::UnitQuaternion{Float64} = UnitQuaternion(1, 0, 0, 0),
		angular_velocity::Vector{Float64} = O,
		Jb = Matrix{Float64}(I, 3, 3),
		kwargs...
		)
	rb = RigidBody(;Jb, invJb=inv(Jb), kwargs...)
	register!(system, rb)

	# Flatten the initial state into the global state vector
	qv = [orientation.w, orientation.x, orientation.y, orientation.z]
	rb_state = [position; velocity; qv; angular_velocity]
	set_state!(system, rb, rb_state)

	return rb
end

function get_position(system::System, rb::RigidBody;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	return rb_state[1:3]
end

function get_velocity(system::System, rb::RigidBody;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	return rb_state[4:6]
end

function get_orientation(system::System, rb::RigidBody;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	return UnitQuaternion(rb_state[7:10])
end

function get_angular_velocity(system::System, rb::RigidBody;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	return rb_state[11:13]
end

function set_position!(system::System, rb::RigidBody, position::Vector{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	rb_state[1:3] = position
	set_state!(system, rb, rb_state, system_state = system_state)
end

function set_velocity!(system::System, rb::RigidBody, velocity::Vector{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	rb_state[4:6] = velocity
	set_state!(system, rb, rb_state, system_state = system_state)
end

function set_orientation!(system::System, rb::RigidBody,
		orientation::UnitQuaternion{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	qv = [orientation.w, orientation.x, orientation.y, orientation.z]
	rb_state[7:10] = qv
	set_state!(system, rb, rb_state, system_state = system_state)
end

function set_angular_velocity!(system::System, rb::RigidBody,
		angular_velocity::Vector{Float64};
		system_state::Union{Vector{Float64},Nothing} = nothing)
	rb_state = get_state(system, rb, system_state = system_state)
	rb_state[11:13] = angular_velocity
	set_state!(system, rb, rb_state, system_state = system_state)
end

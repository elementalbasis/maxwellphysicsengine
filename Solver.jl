using DifferentialEquations
using OrdinaryDiffEqSymplecticRK
#using LinearAlgebra
#using Printf
#using StaticArrays
#using StatsFuns
#using UUIDs
#import Base: @kwdef



###############################################################################
# Solver
###############################################################################

#=
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

	system.evaluation_counter += 1

	return system_state_flow
end

function ODE_function(system_state::Vector{Float64}, system::System, t::Float64)
	return get_state_flow(system, system_state = system_state)
end

function ODE_verlet_function(a_state, v_state, q_state, system::System, t)
    reconstructed_state = reconstruct_state(system, q_state, v_state)

    ia = 1

    for entity in system.entities
        haskey(system.index_map, entity) || continue

        a_entity = get_a_state(system, entity,
                               system_state = reconstructed_state)

        na = length(a_entity)

        if na > 0
            a_state[ia:ia+na-1] .= a_entity
            ia += na
        end
    end

    system.evaluation_counter += 1

    return nothing
end
=#

function ode_function!(a::Vector{Float64}, v::Vector{Float64}, q::Vector{Float64}, system::System, t::Float64)
	state = State(q = q, v = v, n = length(q))

	for entity in system.entities
		if entity isa Thermometer
			update_temperature!(system, entity, state = state)
		end
		degrees_of_freedom(entity) == 0 && continue
		no_qstate(entity) && reset_qstate!(system, entity, state = state)

		range = system.index_map[entity]
		#a[range] .= get_acceleration(system, entity, state = state)
		a[range] = get_astate(system, entity, state = state)
	end

	system.evaluation_counter += 1
end

function update!(system::System, time_increment)
	tspan = (system.time, system.time + time_increment)
	#q_state = get_q_state(system)
	#v_state = get_v_state(system)
	q = system.state.q
	v = system.state.v
	#prob = ODEProblem(ODE_function, system.state, tspan, system)
	#sol = solve(prob, VelocityVerlet(), dt = system.dt)
	#sol = solve(prob)

	prob = SecondOrderODEProblem(ode_function!, v, q, tspan, system)
	sol = solve(prob, VelocityVerlet(), dt = system.dt)

	system.time += time_increment
	#system.state = sol(system.time)
	system.state = get_state_from_solution(system, sol, system.time)

	return sol
end

function get_state_from_solution(system::System, sol, t::Real)
	v = sol(t).x[1]
	q = sol(t).x[2]
	return State(q = q, v = v, n = length(q))
end

#=
function get_state_from_solution(system::System, sol, t::Real)
	v = sol(t).x[1]
	q = sol(t).x[2]
	return reconstruct_state(system, q, v)
end
=#

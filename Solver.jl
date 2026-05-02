using DifferentialEquations
#using LinearAlgebra
#using Printf
#using StaticArrays
#using StatsFuns
#using UUIDs
#import Base: @kwdef



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

function update!(system::System, time_increment)
	tspan = (system.time, system.time + time_increment)
	q_state = get_q_state(system)
	v_state = get_v_state(system)
	#prob = ODEProblem(ODE_function, system.state, tspan, system)
	#sol = solve(prob, VelocityVerlet(), dt = system.dt)
	#sol = solve(prob)
	

	prob = SecondOrderODEProblem(ODE_verlet_function, get_v_state(system), get_q_state(system), tspan, system)
	sol = solve(prob, VelocityVerlet(), dt = system.dt)

	system.time += time_increment
	system.state = sol(system.time)

	return sol
end

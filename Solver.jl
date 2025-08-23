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

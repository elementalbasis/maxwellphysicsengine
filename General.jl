#using DifferentialEquations
#using LinearAlgebra
#using Printf
#using StaticArrays
#using StatsFuns
#using UUIDs
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
abstract type Impulse <: Entity end
abstract type Parameter <: Entity end
entity_state_size(entity::Entity) = 0
degrees_of_freedom(entity::Entity) = 0

@kwdef mutable struct State
	q::Vector{Float64} = []
	v::Vector{Float64} = []
	n::Int64 = 0
	#index_map::Dict{Entity, UnitRange{Int}} = Dict()
end

@kwdef mutable struct System
	forces::Vector{Force} = []
	#impulses::Vector{Impulse} = []
	bodies::Vector{Body} = []
	#parameters::Vector{Parameter} = []
	entities::Vector{Entity} = []

	evaluation_counter::Int64 = 0

	state::State = State()
	time::Float64 = 0.0
	dt::Float64 = 1e-3
	index_map::Dict{Entity, UnitRange{Int}} = Dict()
end

#=
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
=#

function register!(system::System, entity::Entity)
	push!(system.entities, entity)
	entity isa Force && push!(system.forces, entity)
	entity isa Body && push!(system.bodies, entity)

	n = degrees_of_freedom(entity)
	n == 0 && return

	index_start = system.state.n + 1
	index_stop = index_start + n - 1
	index_range = index_start : index_stop
	system.index_map[entity] = index_range

	append!(system.state.q, zeros(n))
	append!(system.state.v, zeros(n))
	system.state.n += n
end

function register!(system::System, entities::Entity...)
	for e in entities
		register!(system, e)
	end
end

#=
function get_q_state(system::System, system_state::Union{Vector{Float64},Nothing} = nothing)
	if system_state == nothing
		system_state = system.state
	end

	q_state = Float64[]
	for entity in system.entities
		haskey(system.index_map, entity) || continue

		q_state = [q_state; get_q_state(system, entity, system_state = system_state)]
	end
	return q_state
end

function get_v_state(system::System, system_state::Union{Vector{Float64},Nothing} = nothing)
	if system_state == nothing
		system_state = system.state
	end

	v_state = Float64[]
	for entity in system.entities
		haskey(system.index_map, entity) || continue

		v_state = [v_state; get_v_state(system, entity, system_state = system_state)]
	end
	return v_state
end

function get_a_state(system::System;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	if system_state == nothing
		system_state = system.state
	end

	a_state = Float64[]

	for entity in keys(system.index_map)
		apppend!(a_state, get_a_state(system, entity, system_state = system_state))
	end

	return a_state
end

function reconstruct_state(system::System, q_state, v_state)
    state = Float64[]

    iq = 1
    iv = 1

    for entity in system.entities
        haskey(system.index_map, entity) || continue

        # probe sizes using current system layout
        q_entity = get_q_state(system, entity)
        v_entity = get_v_state(system, entity)

        nq = length(q_entity)
        nv = length(v_entity)

        if nq > 0
            append!(state, q_state[iq:iq+nq-1])
            iq += nq
        end

        if nv > 0
            append!(state, v_state[iv:iv+nv-1])
            iv += nv
        end
    end

    return state
end
=#

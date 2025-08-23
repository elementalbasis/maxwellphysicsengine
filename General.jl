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
entity_state_size(entity::Entity) = 0

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

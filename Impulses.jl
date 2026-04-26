#using DifferentialEquations
using LinearAlgebra
#using Printf
#using StaticArrays
using StatsFuns
#using UUIDs
import Base: @kwdef



###############################################################################
# Impulses
###############################################################################



# NonOverlap Impulse

@kwdef struct NonOverlap <: Impulse
end

function compute_impulse(system::System, impulse::NonOverlap, particle::Particle;
		system_state::Union{Vector{Float64},Nothing} = nothing)
	return O
end

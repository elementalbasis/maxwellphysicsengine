using DifferentialEquations
using LinearAlgebra
using Printf
using StaticArrays
using StatsFuns
using UUIDs
import Base: @kwdef

include("./General.jl")
include("./Bodies.jl")
include("./Forces.jl")
include("./Solver.jl")

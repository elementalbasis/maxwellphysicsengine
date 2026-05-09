import Pkg

packages = [
	"DifferentialEquations",
	"StaticArrays",
	"StatsFuns",
	"Rotations",
	"ArgParse",
	"Distributions",
	"OrdinaryDiffEqSymplecticRK",
	]

for package in packages
	Pkg.add(package)
end

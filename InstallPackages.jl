import Pkg

packages = [
	"DifferentialEquations",
	"StaticArrays",
	"StatsFuns",
	"Rotations",
	"ArgParse",
	"Distributions",
	]

for package in packages
	Pkg.add(package)
end

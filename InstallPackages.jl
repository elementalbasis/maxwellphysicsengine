import Pkg

packages = [
	"DifferentialEquations",
	"StaticArrays",
	"StatsFuns",
	"Rotations",
	"ArgParse",
	]

for package in packages
	Pkg.add(package)
end

import Pkg

packages = [
	"DifferentialEquations",
	"StaticArrays",
	"StatsFuns",
	"Rotations",
	]

for package in packages
	Pkg.add(package)
end

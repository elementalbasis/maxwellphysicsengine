include("Main.jl")

using Distributions
using ArgParse
s = ArgParseSettings()

global k_boltzmann
@add_arg_table s begin
	"--particle-count"
		arg_type = Int64
		default = 3000
	"--chunk-size"
		arg_type = Float64
		default = 0.1
	"--particle-radius"
		arg_type = Float64
		default = 0.01
	"--particle-mass"
		arg_type = Float64
		default = 1e-3
	"--particle-stiffness"
		arg_type = Float64
		default = 30000
	"--simulation-dt"
		arg_type = Float64
		default = 4e-5
	"--average-kinetic-energy-target"
		arg_type = Float64
		default = 5
	"--cube-dimensions"
		arg_type = Float64
		default = 2
end

args = parse_args(s)

N = args["particle-count"]
r = args["particle-radius"]
m = args["particle-mass"]
k = args["particle-stiffness"]
dt = args["simulation-dt"]
#temperature = args["thermostat-temperature"]
E = args["average-kinetic-energy-target"]
temperature = E * 2/3 / k_boltzmann
L = 1/2 * args["cube-dimensions"]

#energy_target = 0.5
#sensitivity = 0.01

S = System(dt = dt)
H = ChunkGrid(chunk_size = args["chunk-size"])

for i = 1:N
	x = rand(Uniform(-L,L))
	y = rand(Uniform(-L,L))
	z = rand(Uniform(-L,L))
	v = (100 + i/N/1000) * X
	create_particle!(S, radius = r, mass = m, velocity = v, position = x*X+y*Y+z*Z)
end

T = Thermostat(temperature = temperature)
A = ModulatedWallContact(n = X, p = -L*X, k = k, thermostat = T)
B = ModulatedWallContact(n = Y, p = -L*Y, k = k, thermostat = T)
C = ModulatedWallContact(n = Z, p = -L*Z, k = k, thermostat = T)
D = ModulatedWallContact(n = -X, p = L*X, k = k, thermostat = T)
E = ModulatedWallContact(n = -Y, p = L*Y, k = k, thermostat = T)
F = ModulatedWallContact(n = -Z, p = L*Z, k = k, thermostat = T)
G = ParticleContact(k = k, chunk_grid = H)
register!(S, T, A, B, C, D, E, F, G)


T = 120
fps = 7200
#sol = update!(S, T)
delta_t = 0.01

#println("id\tt\tx\ty\tz\tvx\tvy\tvz")
println("t\tE\tEx\tEy\tEz")
while S.time <= T
	sol = update!(S, delta_t)

	for t in range(S.time - delta_t; step = 1/fps, length = Int(delta_t * fps))
		state = get_state_from_solution(S, sol, t)
		Ex = 0
		Ey = 0
		Ez = 0
		E = 0
		for p in S.bodies
			#state = sol[t]
			#x, y, z = get_position(S, p, state = state)
			vx, vy, vz = get_velocity(S, p, state = state)
			Ex += 1/2 * m * vx^2
			Ey += 1/2 * m * vy^2
			Ez += 1/2 * m * vz^2
			#id = p.uuid
			#println(join([id, t, x, y, z, vx, vy, vz], "\t"))
		end

		E = Ex + Ey + Ez

		println(join([t, E, Ex, Ey, Ez], "\t"))
		ux = 3*Ex / E - 1
		uy = 3*Ey / E - 1
		uz = 3*Ez / E - 1
		thr = 0.05
		if (abs(ux) < thr) && (abs(uy) < thr) && (abs(uz) < thr)
			exit()
		end
	end
end

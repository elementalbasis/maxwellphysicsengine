include("Main.jl")

using Distributions

global k_boltzmann
N = 3000
r = 1e-2
m = 1e-3
k = 30000
dt = 4e-5
temperature = 1/k_boltzmann
L = 1

#energy_target = 0.5
#sensitivity = 0.01

S = System(dt = dt)
H = ChunkGrid(chunk_size = 10*r)

for i = 1:N
	#if i == N
	if false
		create_particle!(S, radius = r, mass = m, velocity = 50*Y*i/N)
	else
		x = rand(Uniform(-L,L))
		y = rand(Uniform(-L,L))
		z = rand(Uniform(-L,L))
		create_particle!(S, radius = r, mass = m, velocity = (100 + i/N/1000)*X, position = x*X+y*Y+z*Z)
	end
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


T = 3600
fps = 7200
#sol = update!(S, T)
delta_t = 0.01
i = 0

println("id\tt\tx\ty\tz\tvx\tvy\tvz")
while S.time <= T
	sol = update!(S, delta_t)

	#for t in range(S.time-delta_t, stop = S.time, step = 1/fps)
	global i

	if true
		for t in range(S.time - delta_t; step = 1/fps, length = Int(delta_t * fps))
			state = get_state_from_solution(S, sol, t)
			for p in S.bodies
				#state = sol[t]
				x, y, z = get_position(S, p, state = state)
				vx, vy, vz = get_velocity(S, p, state = state)
				id = p.uuid
				println(join([id, t, x, y, z, vx, vy, vz], "\t"))
			end
		end
	end

	i += 1
end

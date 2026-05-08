include("Main.jl")

N = 3000
r = 1e-2
m = 1e-3
k = 10000
dt = 1e-4

energy_target = 0.5
sensitivity = 0.01

S = System(dt = dt)
H = ChunkGrid(chunk_size = 10*r)

for i = 1:N
	if i == N
		create_particle!(S, radius = r, mass = m, velocity = 10*Y*i/N)
	else
		create_particle!(S, radius = r, mass = m, velocity = (5 + i/N/1000)*X + Z * i/N/1000, position = i/N * Z/1000)
	end
end

A = ModulatedWallContact(n = X, p = -X, k = k, energy_target = energy_target, sensitivity = sensitivity)
B = ModulatedWallContact(n = Y, p = -Y, k = k, energy_target = energy_target, sensitivity = sensitivity)
C = ModulatedWallContact(n = Z, p = -Z, k = k, energy_target = energy_target, sensitivity = sensitivity)
D = ModulatedWallContact(n = -X, p = X, k = k, energy_target = energy_target, sensitivity = sensitivity)
E = ModulatedWallContact(n = -Y, p = Y, k = k, energy_target = energy_target, sensitivity = sensitivity)
F = ModulatedWallContact(n = -Z, p = Z, k = k, energy_target = energy_target, sensitivity = sensitivity)
G = ParticleContact(k = k, chunk_grid = ChunkGrid(chunk_size = 10 * r))
register!(S, A, B, C, D, E, F, G)


T = 3600
fps = 2400
#sol = update!(S, T)
delta_t = 1
i = 0

println("id\tt\tx\ty\tz\tvx\tvy\tvz")
while S.time <= T
	sol = update!(S, delta_t)

	#for t in range(S.time-delta_t, stop = S.time, step = 1/fps)
	global i

	if i % 10 == 0
		for t in range(S.time - delta_t; step = 1/fps, length = delta_t * fps)
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

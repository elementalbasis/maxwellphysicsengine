include("Main.jl")

N = 100
r = 1e-2
m = 1e-3
k = 100

S = System()
H = ChunkGrid(chunk_size = 10*r)

for i = 1:N
	if i == N
		create_particle!(S, radius = r, mass = m, velocity = 10*Y*i/N)
	else
		create_particle!(S, radius = r, mass = m, velocity = (5 + i/N)*X + Z * i/N, position = i/N * Z)
	end
end

A = ModulatedWallContact(n = X, p = -X, k = k, chunk_grid = H)
B = ModulatedWallContact(n = Y, p = -Y, k = k, chunk_grid = H)
C = ModulatedWallContact(n = Z, p = -Z, k = k, chunk_grid = H)
D = ModulatedWallContact(n = -X, p = X, k = k, chunk_grid = H)
E = ModulatedWallContact(n = -Y, p = Y, k = k, chunk_grid = H)
F = ModulatedWallContact(n = -Z, p = Z, k = k, chunk_grid = H)
G = ParticleContact(k = k, chunks = Chunks(cell_size = 100 * r))
register!(S, A, B, C, D, E, F, G)


T = 3600
fps = 1200
#sol = update!(S, T)
delta_t = 0.01

println("id\tt\tx\ty\tz\tvx\tvy\tvz")
while S.time <= T
	sol = update!(S, delta_t)

	for t in range(S.time-delta_t, stop = S.time, step = 1/fps)
		state = sol(t)
		for p in S.bodies
			#state = sol[t]
			x, y, z = get_position(S, p, system_state = state)
			vx, vy, vz = get_velocity(S, p, system_state = state)
			id = p.uuid
			println(join([id, t, x, y, z, vx, vy, vz], "\t"))
		end
	end
end

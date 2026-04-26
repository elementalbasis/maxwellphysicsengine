include("Main.jl")

N = 500
r = 0.01
m = 0.001
k = 100

S = System()

for i = 1:N
	create_particle!(S, radius = r, mass = m, velocity = 10*X*i/N, position = Z * i/N - Y*i/N)
end

A = WallContact(n = X, p = -X, k = k)
B = WallContact(n = Y, p = -Y, k = k)
C = WallContact(n = Z, p = -Z, k = k)
D = WallContact(n = -X, p = X, k = k)
E = WallContact(n = -Y, p = Y, k = k)
F = WallContact(n = -Z, p = Z, k = k)
G = ParticleContact(k = k)
register!(S, A, B, C, D, E, F, G)


T = 3600 * 10
fps = 1200
#sol = update!(S, T)
delta_t = 0.1

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

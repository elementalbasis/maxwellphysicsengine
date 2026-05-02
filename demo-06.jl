include("Main.jl")

r = 0.01
m = 1
k = 100000

S = System()

pa = create_particle!(S, radius = r, mass = m, velocity = -2*X, position = X)
pb = create_particle!(S, radius = r, mass = m, velocity = X, position = -X)

F = ParticleContact(k = k)
register!(S, F)



T = 1
fps = 600
#sol = update!(S, T)
delta_t = 1

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

include("Main.jl")


N = 100
r = 0.1
m = 0.001

S = System()
pa = create_particle!(S, radius = 0.1, position = X, velocity = -X)
pb = create_particle!(S, radius = 0.1, position = -X, velocity = X)

G = ParticleContact()
register!(S, G)


T = 3600 * 10
fps = 60
#sol = update!(S, T)
delta_t = 0.1

println("id\tt\tx\ty\tz")
while S.time <= T
	sol = update!(S, delta_t)

	for t in range(S.time-delta_t, stop = S.time, step = 1/fps)
		state = sol(t)
		for p in S.bodies
			#state = sol[t]
			x, y, z = get_position(S, p, system_state = state)
			id = p.uuid
			println(join([id, t, x, y, z], "\t"))
		end
	end
end

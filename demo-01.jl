include("Mektis.jl")

S = System()
c = create_particle!(S, is_stationary = true)
p = create_particle!(S, mass = 1.0)

G = UniformGravity()
D = LinearDrag()
F = ModulatedSpring(targets = (c, p))
register!(S, G, D, F)

T = 3600
fps = 60
delta_t = 1

println("t\tx\tv")

while S.time <= T

	sol = update!(S, delta_t)
	for t in range(start = S.time - delta_t, stop = S.time, step = 1/fps)
		x = dot(Z, get_position(S, p, system_state = sol(t)))
		v = dot(Z, get_velocity(S, p, system_state = sol(t)))

		println(join([t, x, v], "\t"))
	end
end

include("Main.jl")

S = System()
c = create_particle!(S, is_stationary = true)
p = create_particle!(S, mass = 1.0)

G = UniformGravity()
#D = LinearDrag()
F = Spring(targets = (c, p))
#register!(S, G, D, F)
register!(S, F, G)

T = 60
fps = 60
#delta_t = 1

println("t\tz\tv")

#=
while S.time <= T

	sol = update!(S, delta_t)
	for t in range(start = S.time - delta_t, stop = S.time, step = 1/fps)
		system_state = get_state_from_solution(S, sol, t)
		x = dot(Z, get_position(S, p, system_state = system_state))
		v = dot(Z, get_velocity(S, p, system_state = system_state))

		println(join([t, x, v], "\t"))
	end
end
=#

sol = update!(S, T)
for t in range(start = 0, stop = T, step = 1/fps)
    state = get_state_from_solution(S, sol, t)
    z = get_position(S, p, state = state)[3]
    v = get_velocity(S, p, state = state)[3]
    println(join([t, z, v], "\t"))
end

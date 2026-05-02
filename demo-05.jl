include("Main.jl")

using ArgParse
s = ArgParseSettings()

@add_arg_table s begin
	"--particle-count"
		arg_type = Int64
		default = 500
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
		default = 100
end

args = parse_args(s)

N = args["particle-count"]
r = args["particle-radius"]
m = args["particle-mass"]
k = args["particle-stiffness"]

S = System()
H = ChunkGrid(chunk_size = args["chunk-size"])

for i = 1:N
	if i == N
		create_particle!(S, radius = r, mass = m, velocity = 10*Y*i/N)
	else
		create_particle!(S, radius = r, mass = m, velocity = (5 + i/N)*X + Z * i/N, position = i/N * Z)
	end
end

A = ModulatedWallContact(n = X, p = -X, k = k)
B = ModulatedWallContact(n = Y, p = -Y, k = k)
C = ModulatedWallContact(n = Z, p = -Z, k = k)
D = ModulatedWallContact(n = -X, p = X, k = k)
E = ModulatedWallContact(n = -Y, p = Y, k = k)
F = ModulatedWallContact(n = -Z, p = Z, k = k)
G = ParticleContact(k = k, chunk_grid = H)
register!(S, A, B, C, D, E, F, G)



#sol = update!(S, T)
t = 1

#println("id\tt\tx\ty\tz\tvx\tvy\tvz")

sol = update!(S, t)

state = sol(t)
#=
for p in S.bodies
	#state = sol[t]
	x, y, z = get_position(S, p, system_state = state)
	vx, vy, vz = get_velocity(S, p, system_state = state)
	id = p.uuid
	#println(join([id, t, x, y, z, vx, vy, vz], "\t"))
end
=#

include("Mektis.jl")

N_x = 30
N_y = 30
R = 3
beta = 10
k_max = 1000
sensitivity = 1

S = System()

function correct_bounds(x, N)
	y = mod(x, N)
	return y==0 ? N : y
end

function create_top_ring(S::System)
	m_top = []
	for i_x = 1:N_x
		th = 2*pi*i_x/N_x
		p = create_particle!(S, position = R*X*cos(th)+R*Y*sin(th), is_stationary=true)
		m_top = [m_top; p]
	end
	return m_top
end



function create_body(S::System)
	m_body = []
	for i_y = 1:N_y
		current_layer = []
		for i_x = 1:N_x
			p = create_particle!(S)
			current_layer = [current_layer; p]
		end

		if length(m_body) == 0
			m_body = current_layer
		else
			m_body = [m_body current_layer]
		end
	end
	return m_body
end

function create_bottom_ring(S::System)
	m_bottom = []
	for i_x = 1:N_x
		th = 2*pi*i_x/N_x
		p = create_particle!(S, position = R/2*X*cos(th)+R/2*Y*sin(th);
				    motion_lock = Z)
		m_bottom = [m_bottom; p]
	end
	return m_bottom
end

function connect_particles(S::System, m_top, m_body, m_bottom)
	for i_x = 1:N_x
		pa = m_top[i_x]
		pb = m_body[i_x, 1]
		F = ModulatedSpring(targets = (pa, pb), k_max = k_max, sensitivity = sensitivity)
		register!(S, F)
	end
#=
	for i_y = 1:N_y
		for i_x = 1:N_x
			pa = m_body[i_x, i_y]
			pb = m_body[correct_bounds(i_x+1,N_x), i_y]
			F = ModulatedSpring(pa = pa, pb = pb)
			register!(S, F)
		end
	end
=#
	for i_y = 1:(N_y-1)
		for i_x = 1:N_x
			pa = m_body[i_x,i_y]
			pb = m_body[i_x,i_y+1]
			F = ModulatedSpring(targets=(pa, pb), k_max = k_max, sensitivity = sensitivity)
			register!(S, F)

			pa = m_body[i_x,i_y]
			pb = m_body[correct_bounds(i_x+1,N_x),i_y+1]
			F = ModulatedSpring(targets = (pa, pb), k_max = k_max, sensitivity = sensitivity)
			register!(S, F)
		end
	end

	for i_x = 1:N_x
		pa = m_bottom[i_x]
		pb = m_body[i_x, N_y]
		F = ModulatedSpring(targets=(pa, pb), k_max = k_max, sensitivity = sensitivity)
		register!(S, F)
	end
end

m_top = create_top_ring(S)
m_body = create_body(S)
m_bottom = create_bottom_ring(S)
connect_particles(S, m_top, m_body, m_bottom)

G = UniformGravity()
D = LinearDrag(b = beta)
register!(S, G)
register!(S, D)


T = 3600 * 10
fps = 60
#sol = update!(S, T)
delta_t = 0.1

println("t\tx\ty\tz")
while S.time <= T
	sol = update!(S, delta_t)

	for t in range(S.time-delta_t, stop = S.time, step = 1/fps)
		state = sol(t)
		for p in S.bodies
			#state = sol[t]
			x, y, z = get_position(S, p, system_state = state)
			println(join([t, x, y, z], "\t"))
		end
	end
end

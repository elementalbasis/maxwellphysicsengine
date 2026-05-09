
@kwdef mutable struct Thermostat <: Parameter
	temperature = 300
	sensitivity = 1e-6
end
degrees_of_freedom(thermostat::Thermostat) = 1
no_qstate(thermostat::Thermostat) = true

function get_astate(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	#dT = get_temperature_derivative(system, state = state)
	#return [dT / thermostat.temperature]
	T = get_temperature(system, state = state)
	diff = (T - thermostat.temperature) / thermostat.temperature
	return [diff]
end

#=
function initial_qstate(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	return O
end

function initial_vstate(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	T = get_temperature(system, state = state)
	return [(T - thermostat.temperature) / thermostat.temperature]
end

function set_initial_conditions!(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	qstate = initial_qstate(system, thermostat, state = state)
	vstate = initial_vstate(system, thermostat, state = state)
	set_qstate!(system, thermostat, qstate, state = state)
	set_vstate!(system, thermostat, vstate, state = state)
end

function calibrate_vstate!(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	vstate = initial_vstate(system, thermostat, state = state)
	set_vstate!(system, thermostat, vstate, state = state)
end
=#

function get_thermostat_activation(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	return get_vstate(system, thermostat, state = state)[1] * thermostat.sensitivity
end

function get_thermostat_signal(system::System, thermostat::Thermostat;
		state::Union{State,Nothing} = nothing)
	A = get_thermostat_activation(system, thermostat, state = state)
	#return logistic(A / thermostat.sensitivity)
	return logistic(A)
end




@kwdef mutable struct Thermometer <: Entity
	temperature::Float64 = 0.0
	f::Float64 = 3.0
end
degrees_of_freedom(thermometer::Thermometer) = 0

function update_temperature!(system::System, thermometer::Thermometer;
		state::Union{State,Nothing} = nothing)
	E = get_average_kinetic_energy(system, state = state)
	global k_boltzmann
	T = 2 * E / k_boltzmann / thermometer.f
	thermometer.temperature = T
end

function create_thermometer!(system::System; f::Float64 = 3.0)
	thermometer = Thermometer(f = f)
	register!(S, thermometer)
	update_temperature!(system, thermometer)
	return thermometer
end

function get_temperature(system::System, thermometer::Thermometer;
		state::Union{State,Nothing} = nothing)
	return thermometer.temperature
end

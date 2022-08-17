"""
    kfoldperm(N, k; type="random", breaks=[])

Divides a list of indices 1:N in k groups.

# Arguments
- `N`: length of the list
- `k`: number of groups

# Keywords
- `type`: can be "random" or "ordered"
- `breaks`: ???

# Returns
- `y`: ???
"""
function kfoldperm(N, k; type="random", breaks=[])
	n, r = divrem(N, k)
	b = collect(1:n:N+1)
	Nb = length(b)
	for i in 1:Nb
		b[i] += i > r ? r : i-1
	end
	b = sort(unique(append!(b, breaks)))
	Nbreaks = length(b) - Nb
	if type=="random"
		p = randperm(N)
	elseif type=="ordered"
		p = 1:N
	end
	[p[r] for r in [b[i]:b[i+1]-1 for i=1:k+Nbreaks]] #TODO: use RF starts and ends differently to remove PATCH in run_sim_time_iter
end

"""
    t = points_from_key_times(times; dt)

Returns a vector which contains the same points as `times` but with additional points that
have a separation of at least `dt`.

!!! note
    The last time points could not be present in the output in some cases.
    Some time points could be duplicated in the output.
    Duplicated time points should be removed afterwards (done by
        [`get_variable_times`](@ref)).
    The output represents time points that cannot be missed out during the simulation.

# Arguments
- `times`: (::Vector{Float64}) input time array with key points you want to keep

# Keywords
- `dt`: (::Float64) maximum delta time separation between two time samples

# Returns
- `t`: (::Vector{Float64}) output time array with the same points as the input array but
    with additional points that have a separation of at least `dt`.

# Examples
```julia-repl
julia> times = [0 1 2 10 11 12 20]
1×7 Matrix{Int64}:
 0  1  2  10  11  12  20

julia> points_from_key_times(times; dt=0.5)'
1×46 adjoint(::Vector{Float64}) with eltype Float64:
 0.0  0.5  1.0  1.0  1.5  2.0  2.0  2.5  3.0  3.5  4.0  4.5  …  16.0  16.5  17.0  17.5
   18.0  18.5  19.0  19.5  20.0

julia> points_from_key_times(times; dt=0.7)'
1×32 adjoint(::Vector{Float64}) with eltype Float64:
 0.0  0.7  1.0  1.7  2.0  2.7  3.4  4.1  4.8  5.5  6.2  6.9  …  14.1  14.8  15.5  16.2
   16.9  17.6  18.3  19.0  19.7

julia> points_from_key_times(times; dt=2.5)'
1×16 adjoint(::Vector{Float64}) with eltype Float64:
 0.0  1.0  1.0  2.0  2.0  4.5  7.0  9.5  10.0  11.0  11.0  12.0  12.0  14.5  17.0  19.5
```
"""
function points_from_key_times(times; dt)
    # Fill the `t` empty vector in the `for` loop
	t = Float64[]
	for i = 1:length(times)-1
		if dt < times[i+1] - times[i]
			taux = collect(range(times[i], times[i+1]; step=dt))
		else
			taux = [times[i], times[i+1]]
		end
		append!(t, taux)
	end
	return t
end

"""
    (t, Δt) = get_uniform_times(seq, Δt; Δt_rf=1e-4)

This function, despite its name, actually gets non-uniform time points. Refer to
[`get_variable_times`](@ref) for more details.

!!! note
    This function should be deprecated and the simulator should only use the
    [`get_variable_times`](@ref) function. Note that in this version, this function is
    bypassed by [`get_variable_times`](@ref).

# Arguments
- `seq`: (::Sequence) the sequence object
- `Δt`: (::Real) the nominal delta time separation between two time samples for ADC
    acquisition and Gradients (by nominal we mean that the time separation should be at
    least `Δt` when the samples are regarded by [`is_ADC_on`](@ref) or [`is_GR_on`](@ref)),
    otherwise the time points are not necessary and the separation will be bigger)

# Keywords
- `Δt_rf`: (::Real) the nominal delta time separation between two time samples for RF
    excitation (by nominal we mean that the time separation should be at least `Δt_rf` when
    the samples are regarded by [`is_RF_on`](@ref), otherwise the time points are not
    necessary and the separation will be bigger)

# Returns
- `t`: (::Vector{Float64}) the time array with the non-uniform time values
- `Δt`: (::Vector{Float64}) the delta time array with the separation between two adjacent
    time points of the `t` array
"""
function get_uniform_times(seq, Δt; Δt_rf=1e-4)
	t, Δt = get_variable_times(seq; dt=Δt, dt_rf=Δt_rf)
end

"""
    (t, Δt) = get_variable_times(seq; dt=1, dt_rf=1e-4)

This function returns non-uniform time points that are relevant in the sequence `seq`. It is
important to use a variable time step (instead of constant sampling time) to increase the
simulation speed.

# Arguments
- `seq`: (::Sequence) the sequence object
- `dt`: (::Real) the nominal delta time separation between two time samples for ADC
    acquisition and Gradients (by nominal we mean that the time separation should be at
    least `dt` when the samples are regarded by [`is_ADC_on`](@ref) or [`is_GR_on`](@ref)),
    otherwise the time points are not necessary and the separation will be bigger)

# Keywords
- `Δt_rf`: (::Real) the nominal delta time separation between two time samples for RF
    excitation (by nominal we mean that the time separation should be at least `Δt_rf` when
    the samples are regarded by [`is_RF_on`](@ref), otherwise the time points are not
    necessary and the separation will be bigger)

# Returns
- `t`: (::Vector{Float64}) the time array with the non-uniform time values
- `Δt`: (::Vector{Float64}) the delta time array with the separation between two adjacent
    time points of the `t` array
"""
function get_variable_times(seq; dt=1, dt_rf=1e-4)
	t = Float64[]
	ΔT = durs(seq) #Duration of sequence block
	T0 = cumsum([0; ΔT[:]]) #Start time of each block
	for i = 1:length(seq)
		s = seq[i] #Current sequence block
		t0 = T0[i]
		if is_ADC_on(s)
			ts = get_sample_times(s) .+ t0
			taux = points_from_key_times(ts; dt) # ADC sampling
			append!(t, taux)
		end
		if is_RF_on(s)
			y = s.RF[1]
			delay, T = y.delay, y.T
			t1 = t0 + delay
			t2 = t1 + sum(T)
			rf0 = t0 + get_RF_center(y) #get_RF_center includes delays
			taux = points_from_key_times([t1,rf0,t2]; dt=dt_rf) # Arbitrary RF
			append!(t, taux)
		end
		if is_GR_on(s)
			active_gradients = []
			if is_Gx_on(s) append!(active_gradients, s.GR.x) end
			if is_Gy_on(s) append!(active_gradients, s.GR.y) end
			if is_Gz_on(s) append!(active_gradients, s.GR.z) end
			for y = active_gradients
				ts = get_theo_t(y) .+ t0
				taux = points_from_key_times(ts; dt)
				append!(t, taux)
			end
		end
	end
	t = sort(unique(x -> round(x*1e8), t)) #Removing repeated points
	Δt = t[2:end] .- t[1:end-1]
	t = t[1:end-1]
	t, Δt
end

"""
    get_breaks_in_RF_key_points(seq::Sequence, t)

Calculate RF key points (start-end) to split simulation in RF and non-RF parts.

# Arguments
- `seq::Sequence`: the input sequence object
- `t`: the time array

# Returns
- `key_idxs`: array of key points indexes
"""
function get_breaks_in_RF_key_points(seq::Sequence, t)
	T0 = cumsum([0; durs(seq)[:]])
	# Identifying RF key points
	key_points = Float64[]
	key_idxs = Int[]
	for (i, s) = enumerate(seq)
		if is_RF_on(s)
			t0 = T0[i] + s.RF.delay[1]	#start of RF waverform
			tf = T0[i] + s.RF.dur[1]	#end of RF waveform
			append!(key_points, [t0; tf])
			idx0 = argmin(abs.(t.-t0))
			idxf = argmin(abs.(t.-tf))
			append!(key_idxs,   [idx0; idxf])
		end
	end
	key_idxs
end

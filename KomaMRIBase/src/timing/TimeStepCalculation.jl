const MIN_RISE_TIME = 1e-10

"""
    array_of_ranges = kfoldperm(N, k; breaks=[])

Divides a list of indices from 1 to `N` into `k` groups.

# Arguments
- `N`: (`::Integer`) number of elements to be ordered
- `k`: (`::Integer`) number of groups to divide the `N` elements

# Keywords
- `breaks`: (`::Vector{Integer}`, `=[]`) array of indices where predefined breakpoints are
    placed

# Returns
- `array_of_ranges`: (`::Vector{UnitRange{Integer}}`) array containing ranges of different
    groups. The target is `k` groups, but this could increase by adding elements to the
    `breaks` input array
"""
function kfoldperm(N, k; breaks=[])
	k = min(N, k)
	n, r = divrem(N, k) #N >= k, N < k
	b = collect(1:n:N+1)
	Nb = length(b)
	for i in 1:Nb
		b[i] += i > r ? r : i-1
	end
	b = sort(unique(append!(b, breaks)))
	Nbreaks = length(b) - Nb
	p = 1:N
	return [p[r] for r in [b[i]:b[i+1]-1 for i=1:k+Nbreaks]] #TODO: use RF starts and ends differently to remove PATCH in run_sim_time_iter
end

"""
    t = points_from_key_times(times; dt)

Generates a vector that includes the same points as `times`, with the addition of extra
points separated by no more than `dt`.

# Arguments
- `times`: (`::Vector{Real}`, `[s]`) time array with points you want to keep

# Keywords
- `dt`: (`::Real`, `[s]`) maximum time gap between samples

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with the same points as the input but with
    additional points separated by at most `dt`
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
    t, Δt = get_variable_times(seq; Δt=1e-3, Δt_rf=1e-5)

    This function returns non-uniform time points relevant to the sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `Δt`: (`::Real`, `=1e-3`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients (by nominal we mean that the time separation should be
    at most `Δt` when the samples are regarded by [`KomaMRIBase.is_ADC_on`](@ref) or
    [`KomaMRIBase.is_GR_on`](@ref)), otherwise the time points are not necessary and the
    separation will be larger)
- `Δt_rf`: (`::Real`, `=1e-5`, `[s]`) nominal delta time separation between two time
    samples for RF excitation (by nominal we mean that the time separation should be at most
    `Δt_rf` when the samples are regarded by [`KomaMRIBase.is_RF_on`](@ref), otherwise the
    time points are not necessary and the separation will be larger)

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with non-uniform time values
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time array indicating the separation between two
    adjacent time points in the `t` time array
"""
function get_variable_times(seq; Δt=1e-3, Δt_rf=1e-5)
	t = Float64[]
	ϵ = MIN_RISE_TIME #Small Float64
	T0 = get_block_start_times(seq)
	for i = 1:length(seq)
		s = seq[i] #Current sequence block
		t0 = T0[i]
		if is_RF_on(s)
			y = s.RF[1]
			delay, T = y.delay, y.T
			t1 = t0 + delay
			t2 = t1 + sum(T)
			rf0 = t0 + get_RF_center(y) #get_RF_center includes delays
			taux = points_from_key_times([t1,t1+ϵ,rf0,t2-ϵ,t2]; dt=Δt_rf) # Arbitrary RF. Points (t1+ϵ, t2-ϵ) added to fix bug with ADCs
			append!(t, taux)
		end
		if is_GR_on(s)
			active_gradients = []
			if is_Gx_on(s) append!(active_gradients, s.GR.x) end
			if is_Gy_on(s) append!(active_gradients, s.GR.y) end
			if is_Gz_on(s) append!(active_gradients, s.GR.z) end
			for y = active_gradients
				ts = get_theo_t(y) .+ t0
				taux = points_from_key_times([ts[1]+ϵ; ts; ts[end]-ϵ]; dt=Δt) #The ±ϵ fixes #
				append!(t, taux)
			end
		end
	end
	#Adding ADC samples, and removing repeated points
	tadc = get_adc_sampling_times(seq)
	t = sort(unique([t; tadc])) #Removing repeated points
	#Fixes a problem with ADC at the start and end of the seq
	t0 = t[1]   - ϵ
	tf = t[end] + ϵ
	t = [t0; t; tf]
	#Final time points
	Δt = t[2:end] .- t[1:end-1]
	t, Δt
end

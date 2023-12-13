const MIN_RISE_TIME = 1e-10 #100*eps(1.0) #SmallestFloat
"""
    array_of_ranges = kfoldperm(N, k; type="random", breaks=[])

Divides a list of indices 1:`N` (which is in your imagination) into `k` groups.

!!! note
    It is possible to predifine some break points at specific indices with the `breaks`
    keyword, in this case the number of groups could increase. This is useful to define
    start and end indices of RF pulses to separate the simulation into excitation and
    precession computations.

# Arguments
- `N`: (`::Int64`) the number of elements to be ordered (of an imaginary array 1:`N`)
- `k`: (`::Int64`) the number of groups to divide the `N` elements

# Keywords
- `type`: (`::String`, `="random"`, opts: [`"random"`, `"ordered"`]) the order type option.
    If random, then the indices of the groups are unordered. If "ordered", then the indices
    of the groups are sorted in an incremental order
- `breaks`: (`::Vector{Int64}`, `=[]`) the array of indices where predefined break points
    are placed

# Returns
- `array_of_ranges`: (`::Vector{UnitRange{Int64}}`) the array that contains ranges of
    different groups (the aim target are `k` groups, but this could be increased by adding
    elements in the `breaks` input array)

# Examples
``` julia-repl
julia> kfoldperm(20, 3; type="ordered")
3-element Vector{UnitRange{Int64}}:
 1:7
 8:14
 15:20

julia> kfoldperm(20, 3; type="ordered", breaks=[3])
4-element Vector{UnitRange{Int64}}:
 1:2
 3:7
 8:14
 15:20
```
"""
function kfoldperm(N, k; type="ordered", breaks=[])
	k = min(N,k)
	n, r = divrem(N, k) #N >= k, N < k
	b = collect(1:n:N+1)
	Nb = length(b)
	for i in 1:Nb
		b[i] += i > r ? r : i-1
	end
	b = sort(unique(append!(b, breaks)))
	Nbreaks = length(b) - Nb
	if type=="random"
		p = Random.randperm(N)
	elseif type=="ordered"
		p = 1:N
	end
	[p[r] for r in [b[i]:b[i+1]-1 for i=1:k+Nbreaks]] #TODO: use RF starts and ends differently to remove PATCH in run_sim_time_iter
end

"""
    t = points_from_key_times(times; dt)

Returns a vector which contains the same points as `times` but with additional points that
have a separation of at most `dt`.

!!! note
    The last time points could not be present in the output in some cases.
    Some time points could be duplicated in the output.
    Duplicated time points should be removed afterwards (done by
        [`get_variable_times`](@ref)).
    The output represents all time points regarded during the simulation with a "nominal"
    `dt` separation between two samples.

# Arguments
- `times`: (`::Vector{Float64}`, `[s]`) time array with key points you want to keep

# Keywords
- `dt`: (`::Float64`, `[s]`) maximum delta time separation between two time samples

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with the same points as the input array but with
    additional points that have a separation of at most `dt`.
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

This function returns non-uniform time points that are relevant in the sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `Δt`: (`::Real`, `=1e-3`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients (by nominal we mean that the time separation should be
    at most `Δt` when the samples are regarded by [`KomaMRI.is_ADC_on`](@ref) or
    [`KomaMRI.is_GR_on`](@ref)), otherwise the time points are not necessary and the
    separation will be bigger)
- `Δt_rf`: (`::Real`, `=1e-5`, `[s]`) nominal delta time separation between two time
    samples for RF excitation (by nominal we mean that the time separation should be at most
    `Δt_rf` when the samples are regarded by [`KomaMRI.is_RF_on`](@ref), otherwise the time
    points are not necessary and the separation will be bigger)

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with non-uniform time values
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time array with the separation between two
    adjacent time points of the `t` time array
"""
function get_variable_times(seq; Δt=1e-3, Δt_rf=1e-5)
	t = Float64[]
	ϵ = MIN_RISE_TIME #Small Float64
	ΔT = durs(seq) #Duration of sequence block
	T0 = cumsum([0; ΔT[:]]) #Start time of each block
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

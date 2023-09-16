const EPS = 1e-10 #100*eps(1.0) #SmallestFloat
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

# DEPRECATED?
"""
    t, Δt = get_uniform_times(seq, Δt; Δt_rf=1e-4)

This function, despite its name, actually gets non-uniform time points. Refer to
[`get_variable_times`](@ref) for more details.

!!! note
    This function should be deprecated and the simulator should only use the
    [`get_variable_times`](@ref) function. Note that in this KomaMRI version, this function
    is bypassed by [`get_variable_times`](@ref).

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `[s]`) nominal delta time separation between two time samples for ADC
    acquisition and Gradients (by nominal we mean that the time separation should be at
    most `Δt` when the samples are regarded by [`KomaMRI.is_ADC_on`](@ref) or
    [`KomaMRI.is_GR_on`](@ref)), otherwise the time points are not necessary and the
    separation will be bigger)

# Keywords
- `Δt_rf`: (`::Real`, `=1e-4`, `[s]`) nominal delta time separation between two time
    samples for RF excitation (by nominal we mean that the time separation should be at most
    `Δt_rf` when the samples are regarded by [`KomaMRI.is_RF_on`](@ref), otherwise the time
    points are not necessary and the separation will be bigger)

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with non-uniform time values
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time array with the separation between two
    adjacent time points of the `t` time array
"""
function get_uniform_times(seq, Δt; Δt_rf=1e-4)
	t, Δt = get_variable_times(seq; dt=Δt, dt_rf=Δt_rf)
end

# DEPRECATED?
"""
    t, Δt = get_variable_times(seq; dt=1, dt_rf=1e-4)

This function returns non-uniform time points that are relevant in the sequence `seq`.

!!! note
    It is important to use a variable time step (instead of constant sampling time) to
    increase the simulation speed.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `dt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients (by nominal we mean that the time separation should be
    at most `Δt` when the samples are regarded by [`KomaMRI.is_ADC_on`](@ref) or
    [`KomaMRI.is_GR_on`](@ref)), otherwise the time points are not necessary and the
    separation will be bigger)

# Keywords
- `Δt_rf`: (`::Real`, `=1e-4`, `[s]`) nominal delta time separation between two time
    samples for RF excitation (by nominal we mean that the time separation should be at most
    `Δt_rf` when the samples are regarded by [`KomaMRI.is_RF_on`](@ref), otherwise the time
    points are not necessary and the separation will be bigger)

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with non-uniform time values
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time array with the separation between two
    adjacent time points of the `t` time array
"""
function get_variable_times(seq; dt=1e-3, dt_rf=1e-5)
	t = Float64[]
	ϵ = EPS #Small Float64
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
			taux = points_from_key_times([t1,t1+ϵ,rf0,t2-ϵ,t2]; dt=dt_rf) # Arbitrary RF. Points (t1+ϵ, t2-ϵ) added to fix bug with ADCs
			append!(t, taux)
		end
		if is_GR_on(s)
			active_gradients = []
			if is_Gx_on(s) append!(active_gradients, s.GR.x) end
			if is_Gy_on(s) append!(active_gradients, s.GR.y) end
			if is_Gz_on(s) append!(active_gradients, s.GR.z) end
			for y = active_gradients
				ts = get_theo_t(y) .+ t0
				taux = points_from_key_times([ts[1]+ϵ; ts; ts[end]-ϵ]; dt) #The ±ϵ fixes #
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

##################
### DEPRECATED ###
##################
#"""
#    key_idxs = get_breaks_in_RF_key_points(seq::Sequence, t)
#
#Return the indices of the `t` time array where are RF key points from the `seq` sequence.
#Thus, it is possible to split the simulation into excitation and precession computations.
#
#!!! note
#    By `RF key points` we mean all the start and end points where the RF excitation takes
#    place with the [`KomaMRI.is_RF_on`](@ref) function.
#
## Arguments
#- `seq`: (`::Sequence`) Sequence struct
#- `t`: (`::Vector{Int64}`, `[s]`) non-uniform time array
#
## Returns
#- `key_idxs`: (`::Vector{Int64}`) array of indices of the `t` time array where are RF key
#    points
#"""
#function get_breaks_in_RF_key_points(seq::Sequence, t)
#	T0 = cumsum([0; durs(seq)[:]])
#	# Identifying RF key points
#	key_points = Float64[]
#	key_idxs = Int[]
#	for (i, s) = enumerate(seq)
#		if is_RF_on(s)
#			t0 = T0[i] + s.RF.delay[1]	#start of RF waverform
#			tf = T0[i] + s.RF.dur[1]	#end of RF waveform
#			append!(key_points, [t0; tf])
#			idx0 = findall(t .== t0)
#			idxf = findall(t .== tf)
#			append!(key_idxs,   [idx0; idxf])
#		end
#	end
#	return key_idxs
#end

function get_sim_ranges(seqd::DiscreteSequence; Nblocks)
	ranges = UnitRange{Int}[]
	ranges_bool = Bool[]
	start_idx_rf_block = 0
	start_idx_gr_block = 0
	#Split 1:N into Nblocks like kfoldperm
	N = length(seqd.Δt)
	k = min(N, Nblocks)
	n, r = divrem(N, k) #N >= k, N < k
	breaks = collect(1:n:N+1)
	for i in eachindex(breaks)
		breaks[i] += i > r ? r : i-1
	end
	breaks = breaks[2:end-1] #Remove borders,
	#Iterate over B1 values to decide the simulation UnitRanges
	for i in eachindex(seqd.Δt)
		if abs(seqd.B1[i]) > 10EPS #TODO: This is needed as the function ⏢ in get_rfs is not very accurate
			if start_idx_rf_block == 0 #End RF block
				start_idx_rf_block = i
			end
			if start_idx_gr_block > 0 #End of GR block
				push!(ranges, start_idx_gr_block:i-1)
				push!(ranges_bool, false)
				start_idx_gr_block = 0
			end
		else
			if start_idx_gr_block == 0 #Start GR block
				start_idx_gr_block = i
			end
			if start_idx_rf_block > 0 #End of RF block
				push!(ranges, start_idx_rf_block:i-1)
				push!(ranges_bool, true)
				start_idx_rf_block = 0
			end
		end
		#More subdivisions
		if i in breaks
			if start_idx_rf_block > 0 #End of RF block
				if length(start_idx_rf_block:i-1) > 1
					push!(ranges, start_idx_rf_block:i-1)
					push!(ranges_bool, true)
					start_idx_rf_block = i
				end
			end
			if start_idx_gr_block > 0 #End of RF block
				if length(start_idx_gr_block:i-1) > 1
					push!(ranges, start_idx_gr_block:i-1)
					push!(ranges_bool, false)
					start_idx_gr_block = i
				end
			end
		end
	end
	#Finishing the UnitRange's
	if start_idx_rf_block > 0
		push!(ranges, start_idx_rf_block:N)
		push!(ranges_bool, true)
	end
	if start_idx_gr_block > 0
		push!(ranges, start_idx_gr_block:N)
		push!(ranges_bool, false)
	end
	#Output
	return ranges, ranges_bool
end

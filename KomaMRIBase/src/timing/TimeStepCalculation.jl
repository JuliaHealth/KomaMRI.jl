const MIN_RISE_TIME = 1e-14

"""
    array_of_ranges = kfoldperm(N, k; breaks=[])

Divides a list of indices from 1 to `N` into `k` groups.

# Arguments
- `N`: (`::Integer`) number of elements to be ordered
- `k`: (`::Integer`) number of groups to divide the `N` elements.

# Keywords
- `breaks`: (`::Vector{<:Integer}`, `=[]`) array of indices where predefined breakpoints are
    placed.

# Returns
- `array_of_ranges`: (`::Vector{UnitRange{<:Integer}}`) array containing ranges of different
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

Returns a vector which contains the same points as `times` but with additional points that
have a separation of at most `dt`.


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
- `motion`: (`::Union{NoMotion, Motion, MotionList}`, `=NoMotion()`) phantom motion,
	from which it may be necessary to extract key time points relevant for the simulation

# Returns
- `t`: (`::Vector{Float64}`, `[s]`) time array with non-uniform time values
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time array with the separation between two
    adjacent time points of the `t` time array
"""
function get_variable_times(seq; Δt=1e-3, Δt_rf=1e-5, motion=NoMotion())
	t = Float64[]
	ϵ = MIN_RISE_TIME # Small Float64
	T0 = get_block_start_times(seq)
	for i = 1:length(seq)
        t_block = Float64[]
		s = seq[i] # Current sequence block
		t0 = T0[i]
		if is_RF_on(s)
			y = s.RF[1]
			delay, T = y.delay, y.T
			t1 = t0 + delay
			t2 = t1 + sum(T)
			tc = t0 + delay + get_RF_center(y)
			# Past ~128 s, `t1 + ϵ` (and `t2 - ϵ`) round back to `t1`/`t2`
			# bit-exactly, collapsing the rise-/fall-edge step markers
			# that `linear_interpolation` relies on downstream. The
			# `deduplicate_knots!(_; move_knots=true)` call perturbs any
			# such bit-equal duplicates by 1 ULP via `nextfloat`, so the
			# step survives at any absolute time. No effect when ϵ
			# arithmetic already produces distinct floats (small t).
			key = sort([t1, t1 + ϵ, tc, t2 - ϵ, t2])
			Interpolations.deduplicate_knots!(key; move_knots=true)
			taux = points_from_key_times(key; dt=Δt_rf)
            append!(t_block, taux)
		end
		if is_GR_on(s)
			active_gradients = []
			if is_Gx_on(s) append!(active_gradients, s.GR.x) end
			if is_Gy_on(s) append!(active_gradients, s.GR.y) end
			if is_Gz_on(s) append!(active_gradients, s.GR.z) end
			for y = active_gradients
				ts = times(y) .+ t0
				taux = points_from_key_times([ts[1] + ϵ; ts; ts[end] - ϵ]; dt=Δt)
                append!(t_block, taux)
			end
		end
        if is_ADC_on(s)
            append!(t_block, times(s.ADC[1]) .+ t0) # get_adc_sampling_times(seq) uses times(s.ADC[1]) .+ T0[i]
        end
        append!(t, t_block)
	end
	add_key_time_points!(t, motion)
	# Ensure t is never empty and contains the start and end of the sequence
	append!(t, [0.0, dur(seq)])
	# Removing repeated points
	sort!(unique!(t))
	# Fixes a problem with ADC at the start and end of the seq.
	# Scale the bookend pad with `eps(t)` so it's a distinct Float64
	# at any t (a fixed `ϵ = 1e-14` rounds back onto `t` itself once
	# `eps(t)/2 > ϵ`, ≈128 s), while keeping `MIN_RISE_TIME` as a
	# floor so we don't drop to subnormal precision near t = 0.
	ϵ_lo = max(MIN_RISE_TIME, 2 * eps(t[1]))
	ϵ_hi = max(MIN_RISE_TIME, 2 * eps(t[end]))
	t = [t[1] - ϵ_lo; t; t[end] + ϵ_hi]
	# Time difference
	Δt = t[2:end] .- t[1:end-1]
	t, Δt
end

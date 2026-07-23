# ===========================================================================
# 3. Sequence Sampling Grid
# ===========================================================================
#
# Build the time support where sequence event waveforms must be evaluated.
# Repeated times are preserved where discontinuities need left/right states.

# -- 3.1. Merge candidate sampling times ------------------------------------
count_sorted_equal(times::Union{AbstractVector,AbstractRange}, t) =
    searchsortedlast(times, t) - searchsortedfirst(times, t) + 1
count_sorted_equal(times, t) = count(==(t), times)

function merge_sampling_times(time_sets...)::Vector{Float64}
    all_times = Float64[]
    sizehint!(all_times, sum(length, time_sets))
    foreach(times -> append!(all_times, times), time_sets)
    sort!(all_times)

    merged_times = similar(all_times, 0)
    i = firstindex(all_times)
    while i <= lastindex(all_times)
        t = all_times[i]
        j = i + 1
        while j <= lastindex(all_times) && all_times[j] == t
            j += 1
        end
        n = 1
        if j - i != 1
            n = 0
            for times in time_sets
                m = count_sorted_equal(times, t)
                m > n && (n = m)
            end
        end
        for _ in 1:n
            push!(merged_times, t)
        end
        i = j
    end
    return merged_times
end

# -- 3.2. Preserve event edges and discontinuities ---------------------------
function event_boundary_sampling_times(times)
    isempty(times) && return collect(times)
    issorted(times) || return event_boundary_sampling_times(sort!(collect(times)))

    event_times = eltype(times)[]
    sizehint!(event_times, 4)
    i = firstindex(times)
    while i <= lastindex(times)
        t = times[i]
        j = i
        while j < lastindex(times) && times[j + 1] == t
            j += 1
        end
        (i == firstindex(times) || j == lastindex(times) || i != j) &&
            append!(event_times, Iterators.repeated(t, j - i + 1))
        i = j + 1
    end
    return event_times
end

# -- 3.3. Collapse redundant waveform storage samples -----------------------
amplitude_roundoff_tol(A) = 1000 * eps(float(maximum(abs, A)))

function simplify_waveform_samples(t, A)
    length(t) <= 3 && return t, A
    keep = trues(length(t))
    amplitude_tol = amplitude_roundoff_tol(A)
    for i in 2:(length(t) - 1)
        t[i + 1] == t[i - 1] && continue
        Ai = A[i - 1] + (A[i + 1] - A[i - 1]) * ((t[i] - t[i - 1]) / (t[i + 1] - t[i - 1]))
        keep[i] = !isapprox(A[i], Ai; rtol=0, atol=amplitude_tol)
    end
    count(keep) < 3 && (keep[end - 1] = true)
    return t[keep], A[keep]
end

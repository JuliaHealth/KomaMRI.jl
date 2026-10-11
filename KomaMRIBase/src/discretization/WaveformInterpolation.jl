# ===========================================================================
# 2. Waveform Evaluation
# ===========================================================================
#
# Evaluate sampled RF, gradient, frequency, or ADC gate waveforms on a chosen
# sampling grid.

function linear_interpolate_samples(samples, t; default=zero(eltype(samples.A)), interpolate=true)
    out = similar(samples.A, typeof(default), length(t))
    fill!(out, default)
    _foreach_interpolation_sample(samples, t, interpolate) do indices, lo, hi, w
        value = lo == hi ? samples.A[lo] : samples.A[lo] + (samples.A[hi] - samples.A[lo]) * w
        for k in indices
            out[k] = value
        end
    end
    return out
end

# Share timestamp lookup between direct evaluation and batched array evaluation.
function _foreach_interpolation_sample(f, samples, t, interpolate)
    (isempty(samples.t) || isempty(samples.A)) && return nothing
    last_sample = min(lastindex(samples.t), lastindex(samples.A))
    sample = firstindex(samples.t)
    i = firstindex(t)
    while i <= lastindex(t)
        ti = t[i]
        j = i
        while j < lastindex(t) && t[j + 1] == ti
            j += 1
        end
        while sample <= last_sample && samples.t[sample] < ti
            sample += 1
        end
        if sample <= last_sample && samples.t[sample] == ti
            sample_end = sample
            while sample_end < last_sample && samples.t[sample_end + 1] == ti
                sample_end += 1
            end
            for k in i:j
                l = interpolate ? min(sample + k - i, sample_end) : sample_end - (j - k)
                l >= sample && f(k:k, l, l, 0)
            end
            sample = sample_end + 1
        elseif interpolate && ti >= first(samples.t) && sample <= last_sample
            lo_time, hi_time = samples.t[sample - 1], samples.t[sample]
            w = (ti - lo_time) / (hi_time - lo_time)
            f(i:j, sample - 1, sample, w)
        end
        i = j + 1
    end
    return nothing
end

# Piecewise cubic Hermite interpolation with finite-difference knot slopes.
function cubic_interpolate_samples(samples, t)
    n = length(samples.t)
    lo = clamp.(searchsortedlast.(Ref(samples.t), t), 1, n - 1)
    hi = lo .+ 1
    prev = max.(lo .- 1, 1)
    next = min.(hi .+ 1, n)
    h = samples.t[hi] .- samples.t[lo]
    s = (t .- samples.t[lo]) ./ h
    m_lo = (samples.A[hi] .- samples.A[prev]) ./ (samples.t[hi] .- samples.t[prev])
    m_hi = (samples.A[next] .- samples.A[lo]) ./ (samples.t[next] .- samples.t[lo])
    s2 = s .* s
    s3 = s2 .* s
    return @. (2s3 - 3s2 + 1) * samples.A[lo] +
              (s3 - 2s2 + s) * h * m_lo +
              (-2s3 + 3s2) * samples.A[hi] +
              (s3 - s2) * h * m_hi
end

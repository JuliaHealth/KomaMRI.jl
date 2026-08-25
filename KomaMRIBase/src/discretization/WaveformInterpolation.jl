# ===========================================================================
# 2. Waveform Evaluation
# ===========================================================================
#
# Evaluate sampled RF, gradient, frequency, or ADC gate waveforms on a chosen
# sampling grid.

function linear_interpolate_samples(samples, t; default=zero(eltype(samples.A)), interpolate=true)
    out = similar(samples.A, typeof(default), length(t))
    fill!(out, default)
    (isempty(samples.t) || isempty(samples.A)) && return out
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
                l >= sample && (out[k] = samples.A[l])
            end
            sample = sample_end + 1
        elseif interpolate && ti >= first(samples.t) && sample <= last_sample
            lo_time, hi_time = samples.t[sample - 1], samples.t[sample]
            w = (ti - lo_time) / (hi_time - lo_time)
            lo = samples.A[sample - 1]
            value = lo + (samples.A[sample] - lo) * w
            for k in i:j
                out[k] = value
            end
        end
        i = j + 1
    end
    return out
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

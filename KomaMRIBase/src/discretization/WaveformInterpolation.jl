# ===========================================================================
# 2. Waveform Evaluation
# ===========================================================================
#
# Evaluate sampled RF, gradient, frequency, or ADC gate waveforms on a chosen
# sampling grid.

function linear_interpolate_samples(samples, t; default=zero(eltype(samples.A)), interpolate=true)
    out = Vector{typeof(default)}(undef, length(t))
    isempty(samples.t) && return fill!(out, default)
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
                out[k] = l >= sample ? samples.A[l] : default
            end
            sample = sample_end + 1
        elseif !interpolate || ti < first(samples.t) || sample > last_sample
            out[i:j] .= default
        else
            lo_time, hi_time = samples.t[sample - 1], samples.t[sample]
            w = (ti - lo_time) / (hi_time - lo_time)
            value = samples.A[sample - 1] + (samples.A[sample] - samples.A[sample - 1]) * w
            out[i:j] .= value
        end
        i = j + 1
    end
    return out
end

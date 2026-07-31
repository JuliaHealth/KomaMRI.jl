module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: DiscreteSequence, RF, ampls, freq_times, freqs,
    linear_interpolate_samples, times, with_adc_start_padding

const ReactantRVector = Union{
    Reactant.AnyTracedRVector,
    Reactant.AbstractConcreteArray{<:Number,1},
}

const ReactantRF = RF{<:ReactantRVector,<:Number,<:Number}

const ReactantDiscreteSequence = DiscreteSequence{
    <:AbstractVector,
    <:AbstractVector,
    <:AbstractVector,
    <:ReactantRVector,
}

function _prepend(value, values)
    out = similar(values, length(values) + 1)
    Reactant.allowscalar() do
        out[1] = value
        for i in eachindex(values)
            out[i + 1] = values[i]
        end
    end
    return out
end

function _prepend_first(values)
    out = similar(values, length(values) + 1)
    Reactant.allowscalar() do
        for i in eachindex(values)
            out[i + 1] = values[i]
        end
        out[1] = out[2]
    end
    return out
end

function _copy_like(template, values)
    out = similar(template, eltype(values), size(values))
    Reactant.allowscalar() do
        for i in eachindex(values)
            out[i] = values[i]
        end
    end
    return out
end

function with_adc_start_padding(values::ReactantDiscreteSequence)
    (isempty(values.ADC) || !first(values.ADC)) && return values
    return DiscreteSequence(
        _prepend_first(values.Gx),
        _prepend_first(values.Gy),
        _prepend_first(values.Gz),
        _prepend_first(values.B1),
        _prepend_first(values.Δf),
        _prepend_first(values.ψ),
        _prepend(false, values.ADC),
        _prepend(false, values.excitation_bool),
        _prepend_first(values.t),
        _prepend(zero(eltype(values.Δt)), values.Δt),
    )
end

function _shape_rf_samples(rf::ReactantRF)
    A = cis(rf.ϕ) .* rf.A
    length(A) == 1 && (A = A[[1, 1]])
    out = similar(A, eltype(A), length(A) + 2)
    Reactant.allowscalar() do
        out[1] = zero(eltype(A))
        for i in eachindex(A)
            out[i + 1] = A[i]
        end
        out[end] = zero(eltype(A))
    end
    return out
end

ampls(rf::ReactantRF; freq_in_phase=false) = begin
    freq_in_phase && !iszero(rf.Δf) &&
        throw(ArgumentError("Reactant RF discretization does not yet support frequency offsets folded into phase."))
    _shape_rf_samples(rf)
end

function times(rf::ReactantRF)
    t = length(rf.A) == 1 ?
        [rf.delay, rf.delay, KomaMRIBase.dur(rf), KomaMRIBase.dur(rf)] :
        collect(range(rf.delay, KomaMRIBase.dur(rf); length=length(rf.A)))
    return length(rf.A) == 1 ? t : [first(t); t; last(t)]
end

freq_times(rf::ReactantRF) =
    [rf.delay, rf.delay, KomaMRIBase.dur(rf), KomaMRIBase.dur(rf)]

function freqs(rf::ReactantRF)
    out = fill(rf.Δf, length(freq_times(rf)))
    out[begin] = zero(rf.Δf)
    out[end] = zero(rf.Δf)
    return out
end

function linear_interpolate_samples(
    samples::NamedTuple{(:t,:A),Tuple{TT,AT}},
    t;
    default=zero(eltype(samples.A)),
    interpolate=true,
) where {TT,AT<:ReactantRVector}
    coefficients = zeros(Float64, length(t), length(samples.A))
    defaults = fill(default, length(t))
    if isempty(samples.t)
        return _copy_like(samples.A, coefficients) * samples.A .+
               _copy_like(samples.A, defaults)
    end

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
                if l >= sample
                    coefficients[k, l] = 1.0
                    defaults[k] = zero(default)
                end
            end
            sample = sample_end + 1
        elseif interpolate && ti >= first(samples.t) && sample <= last_sample
            lo_time, hi_time = samples.t[sample - 1], samples.t[sample]
            weight = (ti - lo_time) / (hi_time - lo_time)
            coefficients[i:j, sample - 1] .= 1 - weight
            coefficients[i:j, sample] .= weight
            defaults[i:j] .= zero(default)
        end
        i = j + 1
    end
    return _copy_like(samples.A, coefficients) * samples.A .+
           _copy_like(samples.A, defaults)
end

end

module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: DiscreteSequence, RF, Sequence, ampls, freq_times, freqs,
    linear_interpolate_samples, sample_sequence, times

const ReactantRVector = Union{
    Reactant.AnyTracedRVector,
    Reactant.AbstractConcreteArray{<:Number,1},
}

const ReactantRF = RF{<:ReactantRVector,<:Number,<:Number}

const ReactantSequence = Sequence{
    GT,
    RT,
    AT,
    DT,
    XT,
    DF,
} where {
    GT,
    RT<:AbstractMatrix{<:ReactantRF},
    AT,
    DT,
    XT,
    DF,
}

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

function _with_adc_start_padding(values::ReactantDiscreteSequence)
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

function sample_sequence(
    seq::ReactantSequence;
    motion=KomaMRIBase.NoMotion(),
    sampling_rule=KomaMRIBase.MaxStepSizeRule(1e-3, 5e-5),
    freq_in_phase=false,
)
    length(seq) == 1 || throw(ArgumentError(
        "Reactant RF sequence discretization currently supports one-block sequences only.",
    ))
    T0 = KomaMRIBase.get_block_start_times(seq)
    global_event_times = KomaMRIBase.merge_sampling_times(
        KomaMRIBase.sequence_boundary_sampling_times(seq),
        KomaMRIBase.motion_sampling_times(seq, motion),
    )
    values = KomaMRIBase.sample_sequence_block(
        seq,
        1;
        sampling_rule,
        motion_times=KomaMRIBase.block_global_event_times(T0, 1, global_event_times),
        freq_in_phase,
    )
    return _with_adc_start_padding(values)
end

function _shape_rf_samples(rf::ReactantRF)
    # Traced amplitudes cannot control output shape, so keep the RF sampling grid fixed even
    # when every amplitude evaluates to zero.
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
    defaults = fill(default, length(t))
    (isempty(samples.t) || isempty(samples.A)) &&
        return _copy_like(samples.A, defaults)

    weight_type = typeof(real(default))
    lo_indices = fill(firstindex(samples.A), length(t))
    hi_indices = copy(lo_indices)
    lo_weights = zeros(weight_type, length(t))
    hi_weights = zeros(weight_type, length(t))
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
                    lo_indices[k] = l
                    hi_indices[k] = l
                    lo_weights[k] = one(weight_type)
                    defaults[k] = zero(default)
                end
            end
            sample = sample_end + 1
        elseif interpolate && ti >= first(samples.t) && sample <= last_sample
            lo_time, hi_time = samples.t[sample - 1], samples.t[sample]
            weight = weight_type((ti - lo_time) / (hi_time - lo_time))
            lo_indices[i:j] .= sample - 1
            hi_indices[i:j] .= sample
            lo_weights[i:j] .= one(weight_type) - weight
            hi_weights[i:j] .= weight
            defaults[i:j] .= zero(default)
        end
        i = j + 1
    end
    return samples.A[lo_indices] .* _copy_like(samples.A, lo_weights) .+
           samples.A[hi_indices] .* _copy_like(samples.A, hi_weights) .+
           _copy_like(samples.A, defaults)
end

end

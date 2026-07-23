# ===========================================================================
# 4.1. Max-Step Sampling Rule
# ===========================================================================

# -- 4.1.1. Rule definition --------------------------------------------------
"""
    MaxStepSizeRule(Δt, Δt_rf; preserve_samples=(:gradients,))

Select sequence sampling times so gradient and RF intervals do not exceed `Δt` and `Δt_rf`.
"""
struct MaxStepSizeRule{T,TRF,TPreserved} <: SamplingRule
    Δt::T
    Δt_rf::TRF
    preserve_samples::TPreserved
end

const VALID_PRESERVED_SAMPLES = (:rf, :gradients)

normalize_preserved_samples(sample::Symbol) = normalize_preserved_samples((sample,))

function normalize_preserved_samples(samples)
    samples = Tuple(samples)
    all(sample -> sample in VALID_PRESERVED_SAMPLES, samples) ||
        throw(ArgumentError("preserve_samples must contain only :rf and/or :gradients."))
    return samples
end

MaxStepSizeRule(Δt, Δt_rf; preserve_samples=(:gradients,)) =
    MaxStepSizeRule(Δt, Δt_rf, normalize_preserved_samples(preserve_samples))

preserved_samples(rule::MaxStepSizeRule) = rule.preserve_samples

# -- 4.1.2. Find sampling times that enforce a maximum interval --------------
const MAX_STEP_TIME_SNAP_TOL = MIN_RISE_TIME / 10

function snap_sampling_time(t, times)
    i = searchsortedfirst(times, t)
    i <= lastindex(times) && abs(times[i] - t) <= MAX_STEP_TIME_SNAP_TOL && return times[i]
    i > firstindex(times) && abs(times[i - 1] - t) <= MAX_STEP_TIME_SNAP_TOL && return times[i - 1]
    return t
end

function unique_added_step_times!(added_times, event_times)
    isempty(added_times) && return added_times
    for i in eachindex(added_times)
        added_times[i] = snap_sampling_time(added_times[i], event_times)
    end
    sort!(added_times)
    i = firstindex(added_times)
    for j in (i + 1):lastindex(added_times)
        if added_times[j] - added_times[i] > MAX_STEP_TIME_SNAP_TOL
            i += 1
            added_times[i] = added_times[j]
        end
    end
    resize!(added_times, i)
    return added_times
end

function append_max_step_sampling_times!(added_times, event_times, waveform_times, maxdt; anchors=())
    (isempty(waveform_times) || !isfinite(maxdt) || maxdt <= 0) && return added_times
    checkpoints = snap_sampling_time.((first(waveform_times), anchors..., last(waveform_times)), Ref(event_times))
    append!(added_times, checkpoints)
    step_times = merge_sampling_times(event_times, checkpoints)
    lo = searchsortedfirst(step_times, first(waveform_times))
    hi = searchsortedlast(step_times, last(waveform_times))
    lo < hi || return added_times
    for j in lo:(hi - 1)
        gap = step_times[j + 1] - step_times[j]
        gap > 0 || continue
        step_count = gap / maxdt
        n = round(Int, step_count)
        isapprox(step_count, n; rtol=0, atol=100eps(max(abs(step_count), 1.0))) || (n = ceil(Int, step_count))
        for k in 1:(n - 1)
            push!(added_times, snap_sampling_time(step_times[j] + k * gap / n, event_times))
        end
    end
    return added_times
end

needs_max_step_sampling_times(times, maxdt) =
    !isempty(times) && isfinite(maxdt) && maxdt > 0 && last(times) - first(times) > maxdt

# -- 4.1.3. Add max-step sampling times for event waveforms -------------------
function additional_sampling_times(rule::MaxStepSizeRule, event_times, context::BlockSamplingContext)
    waveforms = context.waveforms
    added_times = eltype(event_times)[]
    if !isempty(waveforms.rf.t)
        needs_max_step_sampling_times(waveforms.rf.t, rule.Δt_rf) &&
            append_max_step_sampling_times!(added_times, event_times, waveforms.rf.t, rule.Δt_rf; anchors=(context.rf_center,))
        waveforms.Δf.t != waveforms.rf.t && needs_max_step_sampling_times(waveforms.Δf.t, rule.Δt_rf) &&
            append_max_step_sampling_times!(added_times, event_times, waveforms.Δf.t, rule.Δt_rf; anchors=(context.rf_center,))
    end
    for gr_waveform in (waveforms.gx, waveforms.gy, waveforms.gz)
        needs_max_step_sampling_times(gr_waveform.t, rule.Δt) &&
            append_max_step_sampling_times!(added_times, event_times, gr_waveform.t, rule.Δt)
    end
    return unique_added_step_times!(added_times, event_times)
end

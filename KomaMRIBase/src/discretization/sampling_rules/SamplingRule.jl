# ===========================================================================
# 4. Sampling Rule Interface
# ===========================================================================
#
# A sampling rule receives physical event times for one sequence block and
# returns additional times to insert. Rules dispatch on their own type, so
# custom policies can be added without branching in discretize.

# -- 4.1. Rule type ----------------------------------------------------------
abstract type SamplingRule end

# -- 4.2. Block context passed to rules --------------------------------------
struct BlockSamplingContext{TDur,TRF,TGrad,TADC,TWaveforms,TCenter}
    duration::TDur
    rf::TRF
    gradients::TGrad
    adc::TADC
    waveforms::TWaveforms
    rf_center::TCenter
end

# -- 4.3. Rule extension point -----------------------------------------------
"""
    times = additional_sampling_times(rule, event_times, context)

Return extra sampling times selected by `rule`.

New sampling rules subtype `SamplingRule` and implement this method. The
returned times are merged with `event_times`; rules add samples, they do
not remove event boundaries. `context` contains the block duration, RF event,
gradient events, ADC event, event waveforms, and RF center time.
"""
additional_sampling_times(::Nothing, event_times, ::BlockSamplingContext) = similar(event_times, 0)

function additional_sampling_times(rule::SamplingRule, event_times, context::BlockSamplingContext)
    throw(ArgumentError("Sampling rule must implement additional_sampling_times(rule, event_times, context)."))
end

preserved_samples(::SamplingRule) = (:gradients,)
preserved_samples(::Nothing) = (:gradients,)
preserve_sample(rule, sample::Symbol) = sample in preserved_samples(rule)

sampling_time_vector(times::AbstractVector) = times
sampling_time_vector(times) = collect(times)

function select_eval_times(rule, event_times, context::BlockSamplingContext)
    new_times = additional_sampling_times(rule, event_times, context)
    isempty(new_times) && return event_times
    new_times = sampling_time_vector(new_times)
    sort!(new_times)
    return merge_sampling_times(event_times, unique!(new_times))
end

select_sampling_times(rule, event_times, context::BlockSamplingContext) =
    select_eval_times(rule, event_times, context)

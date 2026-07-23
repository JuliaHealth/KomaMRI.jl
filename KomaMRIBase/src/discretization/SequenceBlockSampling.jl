# ===========================================================================
# 5. Sequence Block Sampling
# ===========================================================================
#
# A sequence block is sampled by collecting its RF/gradient/ADC event
# waveforms, choosing physical event times, applying the sampling rule, and
# evaluating the event waveforms into a DiscreteSequence segment.

# -- 5.1. Collect block event waveforms -------------------------------------
function block_event_waveforms(rf::RF, gx::Grad, gy::Grad, gz::Grad, adc::ADC; freq_in_phase=false)
    return (;
        rf=event_samples(rf; freq_in_phase),
        Δf=event_samples(rf, Val(:Δf)),
        gx=event_samples(gx),
        gy=event_samples(gy),
        gz=event_samples(gz),
        adc=event_samples(adc),
    )
end

block_event_waveforms(seq, block; freq_in_phase=false) =
    block_event_waveforms(
        seq.RF[1, block],
        seq.GR[1, block],
        seq.GR[2, block],
        seq.GR[3, block],
        seq.ADC[block];
        freq_in_phase,
    )

# -- 5.2. Choose block event times -------------------------------------------
rf_event_times(rule, waveforms) =
    preserve_sample(rule, :rf) ? waveforms.rf.t : event_boundary_sampling_times(waveforms.rf.t)

function gradient_event_times(block_duration, waveforms, gradients)
    return merge_sampling_times(waveforms.gx.t, waveforms.gy.t, waveforms.gz.t, waveforms.adc.t)
end

function gradient_event_times(rule, times)
    preserve_sample(rule, :gradients) ? times : event_boundary_sampling_times(times)
end

function gradient_event_times(block_duration, waveforms, gradients, rule)
    preserve_sample(rule, :gradients) && return gradient_event_times(block_duration, waveforms, gradients)
    return merge_sampling_times(
        gradient_event_times(rule, waveforms.gx.t),
        gradient_event_times(rule, waveforms.gy.t),
        gradient_event_times(rule, waveforms.gz.t),
        waveforms.adc.t,
    )
end

function block_event_times(block_duration, rf::RF, gradients, waveforms, rule)
    if isempty(waveforms.rf.t)
        return gradient_event_times(block_duration, waveforms, gradients, rule)
    end
    return merge_sampling_times(
        rf_event_times(rule, waveforms),
        event_boundary_sampling_times(waveforms.Δf.t),
        gradient_event_times(rule, waveforms.gx.t),
        gradient_event_times(rule, waveforms.gy.t),
        gradient_event_times(rule, waveforms.gz.t),
        waveforms.adc.t,
    )
end

function block_sampling_context(block_duration, rf::RF, gradients, adc::ADC, waveforms)
    return BlockSamplingContext(
        block_duration,
        rf,
        gradients,
        adc,
        waveforms,
        rf_center_time(rf),
    )
end

function block_sampling_times(seq, block; sampling_rule=MaxStepSizeRule(1e-3, 5e-5), waveforms=block_event_waveforms(seq, block), motion_times=())
    rf = seq.RF[1, block]
    gradients = (seq.GR[1, block], seq.GR[2, block], seq.GR[3, block])
    context = block_sampling_context(seq.DUR[block], rf, gradients, seq.ADC[block], waveforms)
    event_times = block_event_times(seq.DUR[block], rf, gradients, waveforms, sampling_rule)
    isempty(motion_times) || (event_times = merge_sampling_times(event_times, motion_times))
    return select_eval_times(sampling_rule, event_times, context)
end

# -- 5.3. Evaluate event waveforms on the evaluation grid --------------------
function sample_sequence_block(waveforms::NamedTuple, eval_times; rf_center, freq_in_phase=false)
    isempty(eval_times) && return DiscreteSequence(eval_times)
    B1 = linear_interpolate_samples(waveforms.rf, eval_times; default=0.0 + 0.0im)
    Δf = linear_interpolate_samples(waveforms.Δf, eval_times; default=0.0)
    ψ = rf_phase_at_times(waveforms.Δf, eval_times, rf_center)
    if freq_in_phase && !isempty(waveforms.Δf.t)
        B1 .*= cis.(ψ)
        fill!(Δf, 0.0)
        fill!(ψ, 0.0)
    end
    Gx = linear_interpolate_samples(waveforms.gx, eval_times; default=0.0)
    Gy = linear_interpolate_samples(waveforms.gy, eval_times; default=0.0)
    Gz = linear_interpolate_samples(waveforms.gz, eval_times; default=0.0)
    ADCflag = linear_interpolate_samples(waveforms.adc, eval_times; default=false, interpolate=false)
    excitation_bool = isempty(waveforms.rf.t) ? fill(false, length(eval_times) - 1) :
        [eval_times[i] < eval_times[i + 1] && first(waveforms.rf.t) <= eval_times[i] && eval_times[i + 1] <= last(waveforms.rf.t) for i in 1:(length(eval_times) - 1)]
    return DiscreteSequence(Gx, Gy, Gz, B1, Δf, ψ, ADCflag, excitation_bool, eval_times, diff(eval_times))
end

# -- 5.4. Sample one sequence block -----------------------------------------
function sample_sequence_block(block_duration, rf::RF, gx::Grad, gy::Grad, gz::Grad, adc::ADC; sampling_rule=MaxStepSizeRule(1e-3, 5e-5), motion_times=(), freq_in_phase=false)
    waveforms = block_event_waveforms(rf, gx, gy, gz, adc)
    context = block_sampling_context(block_duration, rf, (gx, gy, gz), adc, waveforms)
    event_times = block_event_times(block_duration, rf, (gx, gy, gz), waveforms, sampling_rule)
    isempty(motion_times) || (event_times = merge_sampling_times(event_times, motion_times))
    eval_times = select_eval_times(sampling_rule, event_times, context)
    return sample_sequence_block(waveforms, eval_times; rf_center=context.rf_center, freq_in_phase)
end

function sample_sequence_block(seq::Sequence, block; sampling_rule=MaxStepSizeRule(1e-3, 5e-5), motion_times=(), freq_in_phase=false)
    return sample_sequence_block(
        seq.DUR[block],
        seq.RF[1, block],
        seq.GR[1, block],
        seq.GR[2, block],
        seq.GR[3, block],
        seq.ADC[block];
        sampling_rule,
        motion_times,
        freq_in_phase,
    )
end

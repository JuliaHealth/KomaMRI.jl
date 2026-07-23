"""
    signal = CardiacSignal(; heart_rate=1.0, first_peak=0.0)
    signal = CardiacSignal(; r_peaks)
    signal = CardiacSignal(; rr_intervals, first_peak=0.0)

Cardiac trigger signal defined by a constant heart rate, explicit R-peak times,
or RR intervals. Plain heart-rate values are in hertz; Unitful frequencies such
as `60u"minute^-1"` are also accepted. Times are in seconds.
"""
struct CardiacSignal <: AbstractPhysioSignal
    r_peaks::Vector{Float64}
    period::Union{Nothing,Float64}
end

function CardiacSignal(; heart_rate=nothing, r_peaks=nothing, rr_intervals=nothing, first_peak=0.0)
    sources = count(x -> !isnothing(x), (heart_rate, r_peaks, rr_intervals))
    sources == 1 || error("Specify exactly one of `heart_rate`, `r_peaks`, or `rr_intervals`.")

    first_peak = Float64(to_SI(first_peak))
    isfinite(first_peak) || error("`first_peak` must be finite.")
    if !isnothing(heart_rate)
        period = Float64(_heart_rate_period(heart_rate))
        _validate_rr_intervals((period,))
        return CardiacSignal([first_peak], period)
    end

    if !isnothing(rr_intervals)
        rr_intervals = Float64[to_SI(rr) for rr in rr_intervals]
        _validate_rr_intervals(rr_intervals)
        return CardiacSignal([first_peak; first_peak .+ cumsum(rr_intervals)], nothing)
    end

    r_peaks = Float64[to_SI(peak) for peak in r_peaks]
    isempty(r_peaks) && error("`r_peaks` must contain at least one R peak.")
    all(isfinite, r_peaks) || error("`r_peaks` must be finite.")
    all(>(0), diff(r_peaks)) || error("`r_peaks` must be strictly increasing.")
    return CardiacSignal(r_peaks, nothing)
end

_heart_rate_period(heart_rate) = inv(to_SI(heart_rate))

function _validate_rr_intervals(rr_intervals)
    isempty(rr_intervals) && error("`rr_intervals` must contain at least one interval.")
    all(rr -> isfinite(rr) && rr > 0, rr_intervals) ||
        error("RR intervals must be finite and positive.")
    return nothing
end

function _next_r_peak(signal::CardiacSignal, time)
    if !isnothing(signal.period)
        periods = ceil(Int, (time - first(signal.r_peaks) - PULSEQ_TIME_TOL) / signal.period)
        return max(first(signal.r_peaks) + periods * signal.period, time)
    end

    peak = searchsortedfirst(signal.r_peaks, time - PULSEQ_TIME_TOL)
    peak <= length(signal.r_peaks) || error("CardiacSignal has no R peak at or after $(time) s.")
    return max(signal.r_peaks[peak], time)
end

"""
    physio_seq = resolve_triggers(seq, signal)

Return a copy of `seq` with the physiological wait inserted before each trigger.
The trigger occurs at the next R peak, then its `duration` remains the
delay before the next sequence block. KomaMRI currently resolves only standalone
trigger blocks.
"""
function resolve_triggers(seq::Sequence, signal::CardiacSignal)
    input_channels = unique(
        trigger.channel for extensions in seq.EXT for trigger in extensions if _is_trigger(trigger)
    )
    length(input_channels) <= 1 ||
        error("A single CardiacSignal cannot resolve triggers on multiple channels.")

    physio_seq = copy(seq)
    block_start = 0.0
    for block in eachindex(physio_seq.DUR)
        trigger_index = findfirst(_is_trigger, physio_seq.EXT[block])
        if !isnothing(trigger_index)
            block_seq = physio_seq[block]
            (is_GR_on(block_seq) || is_RF_on(block_seq) || is_ADC_on(block_seq)) &&
                error(
                    "Trigger block $block contains gradient, RF, or ADC events, which scanners reject."
                )

            trigger = physio_seq.EXT[block][trigger_index]
            arm_time = block_start + trigger.delay
            wait = _next_r_peak(signal, arm_time) - arm_time
            physio_seq.EXT[block][trigger_index] = Trigger(
                trigger.type,
                trigger.channel,
                trigger.delay + wait,
                trigger.duration,
            )
            physio_seq.DUR[block] += wait
        end
        block_start += physio_seq.DUR[block]
    end
    return physio_seq
end

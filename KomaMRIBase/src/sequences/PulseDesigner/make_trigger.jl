"""
    seq = build_trigger(channel; kwargs...)

Return a one-block `Sequence` with a Pulseq-style trigger extension. See
`make_trigger` for trigger arguments.

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `seq`: Sequence containing the trigger block.
"""
function build_trigger(channel; sys=Scanner(), kwargs...)
    seq = Sequence(sys)
    addblock!(seq, make_trigger(channel; sys, kwargs...))
    return seq
end

"""
    event = make_trigger(channel; kwargs...)

Return a Pulseq-style input trigger extension.

# Arguments
- `channel`: Trigger channel, `:physio1` or `:physio2`.

# Keywords
- `delay=0.0`: Delay before trigger. [`s`]
- `duration=0.0`: Trigger duration. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `event`: `Trigger` extension event.
"""
function make_trigger(channel; delay=0.0, duration=0.0, sys=Scanner())
    delay = to_SI(delay, SIUnitsDefault())
    duration = to_SI(duration, SIUnitsDefault())
    return trigger_event(trigger_channel(channel), delay, duration, sys)
end

trigger_channel(channel) = error("Trigger channel must be a Symbol, for example :physio1 or :physio2.")
trigger_channel(channel::Symbol) = trigger_channel(Val(channel))
trigger_channel(::Val{:physio1}) = 1
trigger_channel(::Val{:physio2}) = 2
trigger_channel(::Val{channel}) where {channel} =
    error("Unsupported trigger channel `:$channel`; use :physio1 or :physio2.")

function trigger_event(channel, delay, duration, sys)
    duration = duration <= sys.GR_Δt ? sys.GR_Δt : duration
    return Trigger(2, channel, delay, duration)
end

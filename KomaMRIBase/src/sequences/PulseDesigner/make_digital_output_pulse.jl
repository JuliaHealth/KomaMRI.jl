"""
    seq = build_digital_output_pulse(channel; kwargs...)

Return a one-block `Sequence` with a Pulseq-style digital output pulse
extension. See `make_digital_output_pulse` for pulse arguments.

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `seq`: Sequence containing the digital output pulse block.
"""
function build_digital_output_pulse(channel; sys=Scanner(), kwargs...)
    seq = Sequence(sys)
    addblock!(seq, make_digital_output_pulse(channel; sys, kwargs...))
    return seq
end

"""
    event = make_digital_output_pulse(channel; kwargs...)

Return a Pulseq-style digital output pulse extension.

# Arguments
- `channel`: Output channel, `:osc0`, `:osc1`, or `:ext1`.

# Keywords
- `delay=0.0`: Delay before output pulse. [`s`]
- `duration=0.0`: Output pulse duration. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `event`: `Trigger` extension event.
"""
function make_digital_output_pulse(channel; delay=0.0, duration=0.0, sys=Scanner())
    delay = to_SI(delay, SIUnitsDefault())
    duration = to_SI(duration, SIUnitsDefault())
    return digital_output_event(digital_output_channel(channel), delay, duration, sys)
end

digital_output_channel(channel) =
    error("Digital output channel must be a Symbol, for example :osc0, :osc1, or :ext1.")
digital_output_channel(channel::Symbol) = digital_output_channel(Val(channel))
digital_output_channel(::Val{:osc0}) = 1
digital_output_channel(::Val{:osc1}) = 2
digital_output_channel(::Val{:ext1}) = 3
digital_output_channel(::Val{channel}) where {channel} =
    error("Unsupported digital output channel `:$channel`; use :osc0, :osc1, or :ext1.")

function digital_output_event(channel, delay, duration, sys)
    duration = duration <= sys.GR_Δt ? sys.GR_Δt : duration
    return Trigger(1, channel, delay, duration)
end

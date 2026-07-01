"""
    seq = build_delay(delay; sys=Scanner())

Return a one-block delay `Sequence`. See `make_delay` for delay arguments.

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `seq`: Sequence containing the delay block.
"""
function build_delay(delay; sys=Scanner())
    seq = Sequence(sys)
    addblock!(seq, make_delay(delay))
    return seq
end

"""
    delay_event = make_delay(delay)

Return a Pulseq-style delay event.

# Arguments
- `delay`: Delay duration. [`s`]

# Returns
- `delay_event`: Delay event.
"""
function make_delay(delay)
    delay = to_SI(delay, SIUnitsDefault())
    isfinite(delay) && delay >= 0 || error("Delay must be finite and nonnegative.")
    return Delay(delay)
end

"""
    seq = build_block_pulse(flip_angle; kwargs...)

Return a one-block `Sequence` with a Pulseq-style block RF pulse. See
`make_block_pulse` for pulse design keywords.

# Arguments
- `flip_angle`: RF flip angle. [`rad`]

# Returns
- `seq`: Sequence containing the RF block.
"""
function build_block_pulse(flip_angle; sys=Scanner(), kwargs...)
    rf = make_block_pulse(flip_angle; sys, kwargs...)
    seq = Sequence(sys)
    addblock!(seq, rf)
    seq.DUR[end] = ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
    return seq
end

"""
    rf = make_block_pulse(flip_angle; duration, kwargs...)
    rf = make_block_pulse(flip_angle; bandwidth, kwargs...)

Return a Pulseq-style block RF event.

# Arguments
- `flip_angle`: RF flip angle. [`rad`]

# Keywords
- `duration`: RF pulse duration. [`s`]
- `bandwidth=nothing`: RF bandwidth used when `duration` is omitted. [`Hz`]
- `time_bw_product=0.0`: Time-bandwidth product used with `bandwidth`.
- `sys=Scanner()`: Scanner defaults and raster times.
- `freq_offset=0.0`: RF frequency offset. [`Hz`]
- `phase_offset=0.0`: RF phase offset. [`rad`]
- `delay=0.0`: RF delay before RF dead-time adjustment. [`s`]
- `use=Undefined()`: RF use label.

If `duration` is omitted, `bandwidth` sets it to `1 / (4bandwidth)`, or
`time_bw_product / bandwidth` when `time_bw_product > 0`.

# Returns
- `rf`: RF event.
"""
function make_block_pulse(flip_angle; duration=nothing, bandwidth=nothing,
    time_bw_product=0.0, sys=Scanner(), freq_offset=0.0, phase_offset=0.0, delay=0.0,
    use=Undefined())
    flip_angle   = to_SI(flip_angle, SIUnitsDefault())
    duration     = to_SI(duration, SIUnitsDefault())
    bandwidth    = to_SI(bandwidth, SIUnitsDefault())
    freq_offset  = to_SI(freq_offset, SIUnitsDefault())
    phase_offset = to_SI(phase_offset, SIUnitsDefault())
    delay        = to_SI(delay, SIUnitsDefault())
    if duration === nothing
        bandwidth === nothing && error("Supply duration or bandwidth.")
        bandwidth > 0 || error("RF bandwidth must be positive.")
        duration = time_bw_product > 0 ? time_bw_product / bandwidth : 1 / (4bandwidth)
    end
    duration > 0 || error("RF pulse duration must be positive.")
    rf_duration = round_to_raster(duration, sys.RF_Δt)
    rf_duration > 0 || error("RF pulse duration is shorter than the RF raster.")
    amplitude = flip_angle / (2π * γ * duration)
    return RF(amplitude, rf_duration, freq_offset, max(delay, sys.RF_dead_time);
        ϕ=phase_offset, use)
end

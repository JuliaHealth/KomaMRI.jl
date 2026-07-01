"""
    seq = build_arbitrary_rf(signal, flip_angle; kwargs...)

Return a `Sequence` with a Pulseq-style arbitrary RF pulse. See
`make_arbitrary_rf` for pulse design keywords.

# Arguments
- `signal`: RF waveform shape.
- `flip_angle`: RF flip angle. [`rad`]

# Returns
- `seq`: Sequence containing the RF block and optional slice gradients.
"""
function build_arbitrary_rf(signal, flip_angle; sys=Scanner(), kwargs...)
    rf, gz, gz_rephaser = make_arbitrary_rf(signal, flip_angle; sys, kwargs...)
    seq = Sequence(sys)
    if gz === nothing
        addblock!(seq, rf)
        seq.DUR[end] = ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
        return seq
    end
    addblock!(seq, rf; z=gz)
    seq.DUR[end] = ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
    addblock!(seq; z=gz_rephaser)
    return seq
end

"""
    rf, gz, gzr = make_arbitrary_rf(signal, flip_angle; sys=Scanner(), kwargs...)

Return a Pulseq-style arbitrary RF event tuple. `gz` and `gzr` are `nothing`
unless `slice_thickness` is supplied.

# Arguments
- `signal`: RF waveform shape.
- `flip_angle`: RF flip angle. [`rad`]

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.
- `slice_thickness=nothing`: Slice thickness for slice-selective RF. [`m`]
- `bandwidth=nothing`: RF bandwidth for slice selection. [`Hz`]
- `time_bw_product=0.0`: Time-bandwidth product used with `slice_thickness`.
- `freq_offset=0.0`: RF frequency offset. [`Hz`]
- `phase_offset=0.0`: RF phase offset. [`rad`]
- `delay=0.0`: RF delay before RF dead-time adjustment. [`s`]
- `dwell=sys.RF_Δt`: RF sample spacing. [`s`]
- `center=nothing`: RF center relative to the Pulseq shape start. [`s`]
- `use=Undefined()`: RF use label.
- `max_grad=nothing`: Slice-gradient amplitude limit override. Plain numbers use Pulseq units. [`Hz/m`]
- `max_slew=nothing`: Slice-gradient slew limit override. Plain numbers use Pulseq units. [`Hz/m/s`]

# Returns
- `rf`: RF event.
- `gz`: Slice-select gradient event.
- `gzr`: Slice rephaser gradient event.
"""
function make_arbitrary_rf(signal, flip_angle; sys=Scanner(), slice_thickness=nothing,
    bandwidth=nothing, time_bw_product=0.0, freq_offset=0.0, phase_offset=0.0,
    delay=0.0, dwell=sys.RF_Δt, center=nothing, use=Undefined(), max_grad=nothing,
    max_slew=nothing)
    flip_angle      = to_SI(flip_angle, SIUnitsDefault())
    slice_thickness = to_SI(slice_thickness, SIUnitsDefault())
    bandwidth       = to_SI(bandwidth, SIUnitsDefault())
    freq_offset     = to_SI(freq_offset, SIUnitsDefault())
    phase_offset    = to_SI(phase_offset, SIUnitsDefault())
    delay           = to_SI(delay, SIUnitsDefault())
    dwell           = to_SI(dwell, SIUnitsDefault())
    center          = to_SI(center, SIUnitsDefault())
    max_grad        = isnothing(max_grad) ? sys.Gmax : to_SI(max_grad, PulseqUnitsDefault())
    max_slew        = isnothing(max_slew) ? sys.Smax : to_SI(max_slew, PulseqUnitsDefault())
    dwell > 0 || error("RF dwell time must be positive.")
    waveform = vec(collect(signal))
    isempty(waveform) && error("RF signal must not be empty.")
    duration = length(waveform) * dwell
    waveform = normalize_flip_angle(waveform, dwell, flip_angle) .|> complex
    rf_start_time = max(delay, sys.RF_dead_time)
    rf_center_value = center === nothing || !isfinite(center) ? nothing :
        clamp(center, 0, duration) - dwell / 2
    rf = RF(waveform, (length(waveform) - 1) * dwell, freq_offset,
        rf_start_time + dwell / 2; center=rf_center_value, ϕ=phase_offset, use)
    slice_thickness === nothing && return rf, nothing, nothing
    slice_thickness > 0 || error("Slice thickness must be positive.")
    if time_bw_product > 0
        bandwidth === nothing || error("Supply bandwidth or time_bw_product, not both.")
        bandwidth = time_bw_product / duration
    elseif bandwidth === nothing
        error("Supply bandwidth or time_bw_product for slice selection.")
    end
    bandwidth > 0 || error("RF bandwidth must be positive.")
    slice_area = bandwidth * duration / (γ * slice_thickness)
    gz = trapezoid(;
        flat_time=duration, flat_area=slice_area, sys, max_grad, max_slew,
    )
    if rf_start_time > gz.rise
        gz.delay = ceil_to_raster(rf_start_time - gz.rise, sys.GR_Δt)
    end
    if rf_start_time < gz.rise + gz.delay
        rf.delay = gz.rise + gz.delay + dwell / 2
    end
    rephaser_area = slice_rephaser_area(gz, slice_area, rf, duration, sys)
    gz_rephaser = trapezoid(; area=rephaser_area, sys, max_grad, max_slew)
    return rf, gz, gz_rephaser
end

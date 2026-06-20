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
    dwell = length(rf.A) > 1 ? rf.T / (length(rf.A) - 1) : sys.RF_Δt
    rf_end = rf.delay + rf.T + dwell / 2 + sys.RF_ring_down_time
    if gz === nothing
        addblock!(seq, rf, Duration(ceil_to_raster(rf_end, sys.DUR_Δt)))
        return seq
    end
    duration = ceil_to_raster(max(rf_end, dur(gz)), sys.DUR_Δt)
    addblock!(seq, rf, Duration(duration); z=gz)
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
- `max_grad=sys.Gmax`: Slice-gradient amplitude limit override. [`T/m`]
- `max_slew=sys.Smax`: Slice-gradient slew limit override. [`T/m/s`]

# Returns
- `rf`: RF event.
- `gz`: Slice-select gradient event.
- `gzr`: Slice rephaser gradient event.
"""
function make_arbitrary_rf(signal, flip_angle; sys=Scanner(), slice_thickness=nothing,
    bandwidth=nothing, time_bw_product=0.0, freq_offset=0.0, phase_offset=0.0,
    delay=0.0, dwell=sys.RF_Δt, center=nothing, use=Undefined(), max_grad=sys.Gmax,
    max_slew=sys.Smax)
    dwell > 0 || error("RF dwell time must be positive.")
    waveform = vec(collect(signal))
    isempty(waveform) && error("RF signal must not be empty.")
    duration = length(waveform) * dwell
    waveform = normalize_flip_angle(waveform, dwell, flip_angle) .|> complex
    rf_start_time = max(delay, sys.RF_dead_time)
    rf_center = center === nothing || !isfinite(center) ? nothing :
        clamp(center, 0, duration) - dwell / 2
    rf = RF(waveform, (length(waveform) - 1) * dwell, freq_offset,
        rf_start_time + dwell / 2; center=rf_center, ϕ=phase_offset, use)
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
    pulseq_center = rf.center + dwell / 2
    gz = make_trapezoid(;
        flat_time=duration, flat_area=slice_area, sys, max_grad, max_slew,
    )
    if rf_start_time > gz.rise
        gz.delay = ceil_to_raster(rf_start_time - gz.rise, sys.GR_Δt)
    end
    if rf_start_time < gz.rise + gz.delay
        rf.delay = gz.rise + gz.delay + dwell / 2
    end
    # MATLAB Pulseq makeArbitraryRf uses rf.center in seconds in this formula.
    rephaser_area = -slice_area * (1 - pulseq_center) / duration -
        (area(gz) - slice_area) / 2
    gz_rephaser = make_trapezoid(; area=rephaser_area, sys, max_grad, max_slew)
    return rf, gz, gz_rephaser
end

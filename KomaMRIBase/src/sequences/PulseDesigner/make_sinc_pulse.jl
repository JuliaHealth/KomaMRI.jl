"""
    seq = build_sinc_pulse(flip_angle; kwargs...)

Return a `Sequence` with a Pulseq-style sinc RF pulse. See `make_sinc_pulse`
for pulse design keywords.

# Arguments
- `flip_angle`: RF flip angle. [`rad`]

# Returns
- `seq`: Sequence containing the RF block and optional slice gradients.
"""
function build_sinc_pulse(flip_angle; sys=Scanner(), kwargs...)
    rf, gz, gz_rephaser = make_sinc_pulse(flip_angle; sys, kwargs...)
    seq = Sequence(sys)
    if gz === nothing
        addblock!(seq, rf, Duration(block_duration_to_fit_rf_ringdown(rf, sys)))
        return seq
    end
    duration = block_duration_to_fit_rf_ringdown(rf, sys, gz)
    addblock!(seq, rf, Duration(duration); z=gz)
    addblock!(seq; z=gz_rephaser)
    return seq
end

"""
    rf, gz, gzr = make_sinc_pulse(flip_angle; duration, sys=Scanner(), kwargs...)

Return a Pulseq-style sinc RF event tuple. `gz` and `gzr` are `nothing` unless
`slice_thickness` is supplied.

# Arguments
- `flip_angle`: RF flip angle. [`rad`]

# Keywords
- `duration`: RF pulse duration. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.
- `slice_thickness=nothing`: Slice thickness for slice-selective RF. [`m`]
- `freq_offset=0.0`: RF frequency offset. [`Hz`]
- `phase_offset=0.0`: RF phase offset. [`rad`]
- `time_bw_product=4.0`: Time-bandwidth product.
- `apodization=0.0`: Cosine-window weight.
- `center_pos=0.5`: RF center position as a fraction of `duration`.
- `delay=0.0`: RF delay before RF dead-time adjustment. [`s`]
- `dwell=sys.RF_Δt`: RF sample spacing. [`s`]
- `use=Undefined()`: RF use label.
- `max_grad=sys.Gmax`: Slice-gradient amplitude limit override. [`T/m`]
- `max_slew=sys.Smax`: Slice-gradient slew limit override. [`T/m/s`]

# Returns
- `rf`: RF event.
- `gz`: Slice-select gradient event.
- `gzr`: Slice rephaser gradient event.
"""
function make_sinc_pulse(flip_angle; duration, sys=Scanner(), slice_thickness=nothing,
    freq_offset=0.0, phase_offset=0.0, time_bw_product=4.0, apodization=0.0,
    center_pos=0.5, delay=0.0, dwell=sys.RF_Δt, use=Undefined(),
    max_grad=sys.Gmax, max_slew=sys.Smax)
    duration > 0 || error("RF pulse duration must be positive.")
    rf_start_time = max(delay, sys.RF_dead_time)
    rf = sinc_rf_event(
        flip_angle, duration, dwell, time_bw_product, apodization, center_pos;
        rf_start_time, freq_offset, phase_offset, use,
    )
    slice_thickness === nothing && return rf, nothing, nothing
    slice_thickness > 0 || error("Slice thickness must be positive.")
    slice_area = time_bw_product / (γ * slice_thickness)
    gz, gz_rephaser = slice_select_gradient_events(duration, slice_area, center_pos; sys,
        rf_start_time, max_grad, max_slew)
    rf.delay = max(rf_start_time, gz.rise + gz.delay) + dwell / 2
    return rf, gz, gz_rephaser
end

function sinc_rf_event(flip_angle, duration, dwell, time_bw_product, apodization,
    center_pos; rf_start_time, freq_offset, phase_offset, use)
    n = raster_samples(duration, dwell)
    t = ((1:n) .- 0.5) .* dwell
    tt = t .- duration * center_pos
    waveform = @. (1 - apodization + apodization * cos(2π * tt / duration)) *
        sinc(time_bw_product * tt / duration)
    waveform = normalize_flip_angle(waveform, dwell, flip_angle)
    rf = RF(waveform, (n - 1) * dwell, freq_offset, rf_start_time + dwell / 2;
        center=duration * center_pos - dwell / 2, ϕ=phase_offset, use)
    return rf
end

function slice_select_gradient_events(
    duration, slice_area, center_pos; sys, rf_start_time, max_grad, max_slew,
)
    gz = make_trapezoid(;
        flat_time=duration, flat_area=slice_area, sys, max_grad, max_slew,
    )
    rephaser_area = -slice_area * (1 - center_pos) - (area(gz) - slice_area) / 2
    gz_rephaser = make_trapezoid(; area=rephaser_area, sys, max_grad, max_slew)
    if rf_start_time > gz.rise
        gz.delay = ceil_to_raster(rf_start_time - gz.rise, sys.GR_Δt)
    end
    return gz, gz_rephaser
end

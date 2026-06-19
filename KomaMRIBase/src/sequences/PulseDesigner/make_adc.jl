"""
    seq = build_adc(num_samples, dwell; kwargs...)
    seq = build_adc(num_samples; dwell, kwargs...)
    seq = build_adc(num_samples; duration, kwargs...)

Return a one-block ADC `Sequence`. See `make_adc` for ADC arguments.

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `seq`: Sequence containing the ADC block.
"""
function build_adc(num_samples, args...; sys=Scanner(), kwargs...)
    adc = make_adc(num_samples, args...; sys, kwargs...)
    seq = Sequence(sys)
    adc_end = KomaMRIBase._pulseq_duration(adc) + sys.ADC_dead_time
    duration = ceil_to_raster(adc_end, sys.DUR_Δt)
    addblock!(seq, adc, Duration(duration))
    return seq
end

"""
    adc = make_adc(num_samples, dwell; kwargs...)
    adc = make_adc(num_samples; dwell, kwargs...)
    adc = make_adc(num_samples; duration, kwargs...)

Return a Pulseq-style ADC event.

# Arguments
- `num_samples`: Number of ADC samples.

# Keywords
- `dwell=nothing`: ADC dwell time. [`s`]
- `duration=nothing`: ADC sampling duration, equal to `num_samples * dwell`. [`s`]
- `delay=0.0`: Pulseq ADC delay before the first dwell interval. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.
- `freq_offset=0.0`: ADC frequency offset. [`Hz`]
- `phase_offset=0.0`: ADC phase offset. [`rad`]

Exactly one of `dwell` or `duration` must be supplied.

# Returns
- `adc`: ADC event.
"""
function make_adc(num_samples; dwell=nothing, duration=nothing, delay=0.0,
    sys=Scanner(), freq_offset=0.0, phase_offset=0.0)
    return make_adc(num_samples, dwell; duration, delay, sys, freq_offset, phase_offset)
end

function make_adc(num_samples, ::Nothing; duration=nothing, delay=0.0,
    sys=Scanner(), freq_offset=0.0, phase_offset=0.0)
    isinteger(num_samples) || error("num_samples must be integer-valued.")
    num_samples > 0 || error("num_samples must be positive.")
    duration === nothing && error("Supply dwell or duration.")
    duration > 0 || error("ADC duration must be positive.")
    return make_adc(
        num_samples, duration / num_samples; delay, sys, freq_offset, phase_offset,
    )
end

function make_adc(num_samples, dwell; duration=nothing, delay=0.0,
    sys=Scanner(), freq_offset=0.0, phase_offset=0.0)
    isinteger(num_samples) || error("num_samples must be integer-valued.")
    num_samples > 0 || error("num_samples must be positive.")
    duration === nothing || error("Supply dwell or duration, not both.")
    num_samples = Int(num_samples)
    dwell > 0 || error("ADC dwell must be positive.")
    delay = max(delay, sys.ADC_dead_time) + dwell / 2
    sampling_time = num_samples == 1 ? dwell : (num_samples - 1) * dwell
    return ADC(num_samples, sampling_time, delay, freq_offset, phase_offset)
end

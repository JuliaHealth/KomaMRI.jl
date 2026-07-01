"""
    seq = build_extended_trapezoid(axis::Symbol, times, amplitudes; kwargs...)

Return a one-block `Sequence` with a Pulseq-style extended trapezoid on `axis`.
See `make_extended_trapezoid` for gradient design keywords.

# Arguments
- `axis::Symbol`: Gradient axis, one of `:x`, `:y`, or `:z`.
- `times`: Gradient sample times. [`s`]
- `amplitudes`: Gradient amplitudes at `times`. Plain numbers use Pulseq units. [`Hz/m`]

# Returns
- `seq`: Sequence containing the gradient block.
"""
function build_extended_trapezoid(axis::Symbol, times, amplitudes; sys=Scanner(), kwargs...)
    seq = Sequence(sys)
    axis in (:x, :y, :z) || error("Gradient axis must be :x, :y, or :z.")
    grad = make_extended_trapezoid(times, amplitudes; sys, kwargs...)
    addblock!(seq; (; axis => grad)...)
    return seq
end

"""
    grad = make_extended_trapezoid(times, amplitudes; kwargs...)

Return a Pulseq-style extended trapezoid gradient event.

# Arguments
- `times`: Gradient sample times. [`s`]
- `amplitudes`: Gradient amplitudes at `times`. Plain numbers use Pulseq units. [`Hz/m`]

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.
- `max_grad=nothing`: Amplitude limit for `convert2arbitrary=true`. Plain numbers use Pulseq units. [`Hz/m`]
- `max_slew=nothing`: Slew limit for `convert2arbitrary=true`. Plain numbers use Pulseq units. [`Hz/m/s`]
- `skip_check=false`: Skip Pulseq's nonzero-start connection check.
- `convert2arbitrary=false`: If true, sample the points as an arbitrary gradient.

# Returns
- `grad`: Extended trapezoid gradient event.
"""
function make_extended_trapezoid(times, amplitudes; sys=Scanner(), max_grad=nothing,
    max_slew=nothing, skip_check=false, convert2arbitrary=false)
    times      = to_SI.(times, Ref(SIUnitsDefault()))
    amplitudes = to_SI.(amplitudes, Ref(PulseqUnitsDefault()))
    max_grad   = isnothing(max_grad) ? sys.Gmax : to_SI(max_grad, PulseqUnitsDefault())
    max_slew   = isnothing(max_slew) ? sys.Smax : to_SI(max_slew, PulseqUnitsDefault())
    return extended_trapezoid(
        times, amplitudes; sys, max_grad, max_slew, skip_check, convert2arbitrary,
    )
end

function extended_trapezoid(times, amplitudes; sys=Scanner(), max_grad=sys.Gmax,
    max_slew=sys.Smax, skip_check=false, convert2arbitrary=false)
    length(times) == length(amplitudes) ||
        error("Times and amplitudes must have the same length.")
    all(iszero, times) && error("At least one time must be nonzero.")
    any(diff(times) .<= 0) && error("Times must be strictly increasing.")

    rounded_times = round_to_raster.(times, sys.GR_Δt)
    isapprox(times[end], rounded_times[end]; rtol=0, atol=KomaMRIBase.PULSEQ_TIME_TOL) ||
        error("The last time point must be on the gradient raster.")

    if !skip_check
        times[1] > 0 && amplitudes[1] != 0 &&
            error("Nonzero first amplitude must connect to previous block.")
    end

    if convert2arbitrary
        first_tick = round(Int, minimum(times) / sys.GR_Δt)
        last_tick = round(Int, maximum(times) / sys.GR_Δt) - 1
        sample_times = ((first_tick:last_tick) .* sys.GR_Δt) .+ sys.GR_Δt / 2
        waveform = linear_interpolation(times, amplitudes).(sample_times)
        duration = (length(waveform) - 1) * sys.GR_Δt
        delay = first_tick * sys.GR_Δt
        edge_time = sys.GR_Δt / 2
        first_amp = first(amplitudes)
        last_amp = last(amplitudes)
        grad = Grad(waveform, duration, edge_time, edge_time, delay, first_amp, last_amp)
        check_hw_limits(grad; max_grad, max_slew)
        return grad
    end

    all(isapprox.(times, rounded_times; rtol=0, atol=KomaMRIBase.PULSEQ_TIME_TOL)) ||
        error("All time points must be on the gradient raster.")
    delay = rounded_times[1]
    durations = diff(rounded_times .- delay)
    return Grad(amplitudes, durations, 0.0, 0.0, delay, first(amplitudes), last(amplitudes))
end

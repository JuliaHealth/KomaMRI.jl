"""
    seq = build_trapezoid(axis::Symbol; kwargs...)

Return a one-block `Sequence` with a Pulseq-style trapezoidal gradient on
`axis`. See `make_trapezoid` for gradient design keywords.

# Arguments
- `axis::Symbol`: Gradient axis for `build_trapezoid`, one of `:x`, `:y`, or `:z`.

# Returns
- `seq`: Sequence containing the gradient block.
"""
function build_trapezoid(axis::Symbol; sys=Scanner(), kwargs...)
    seq = Sequence(sys)
    axis in (:x, :y, :z) || error("Gradient axis must be :x, :y, or :z.")
    addblock!(seq; (; axis => make_trapezoid(; sys, kwargs...))...)
    return seq
end

"""
    grad = make_trapezoid(; area, kwargs...)
    grad = make_trapezoid(; flat_area, flat_time, kwargs...)
    grad = make_trapezoid(; amplitude, flat_time, kwargs...)
    grad = make_trapezoid(; amplitude, duration, kwargs...)

Return a Pulseq-style trapezoidal gradient event.

# Keywords
- `area=nothing`: Total gradient area, including ramps. [`T/m*s`]
- `flat_area=nothing`: Flat-top area, excluding ramps. [`T/m*s`]
- `amplitude=nothing`: Flat-top gradient amplitude. [`T/m`]
- `flat_time=nothing`: Flat-top duration. [`s`]
- `duration=nothing`: Total gradient duration, excluding `delay`. [`s`]
- `delay=0.0`: Gradient delay. [`s`]
- `rise_time=nothing`: Rise ramp duration. [`s`]
- `fall_time=nothing`: Fall ramp duration; requires `rise_time`. [`s`]
- `sys=Scanner()`: Scanner defaults and raster times.
- `max_grad=sys.Gmax`: Gradient amplitude limit override. [`T/m`]
- `max_slew=sys.Smax`: Slew limit override. [`T/m/s`]

Exactly one of `area`, `flat_area`, or `amplitude` must be supplied.

# Returns
- `grad`: Trapezoidal gradient event.
"""
function make_trapezoid(; area=nothing, flat_area=nothing, amplitude=nothing,
    flat_time=nothing, duration=nothing, delay=0.0, rise_time=nothing, fall_time=nothing,
    sys=Scanner(), max_grad=sys.Gmax, max_slew=sys.Smax)
    # Keyword dispatch cannot distinguish Unitful `area`, `flat_area`, or `amplitude`.
    area = si_gradient_area(area)
    flat_area = si_gradient_area(flat_area)
    amplitude = si_gradient(amplitude)
    flat_time = si_time(flat_time)
    duration = si_time(duration)
    delay = si_time(delay)
    rise_time = si_time(rise_time)
    fall_time = si_time(fall_time)
    max_grad = si_gradient(max_grad)
    max_slew = si_slew_rate(max_slew)

    if count(!isnothing, (area, flat_area, amplitude)) != 1
        error("Supply exactly one of area, flat_area, or amplitude.")
    end
    if area !== nothing
        if duration === nothing
            return shortest_trapezoid(area, delay, sys.GR_Δt, max_grad, max_slew)
        end
        return trapezoid_from_area(area, duration; delay, sys, max_grad, max_slew,
            rise_time, fall_time)
    elseif flat_area !== nothing
        flat_time === nothing && error("Supply flat_time with flat_area.")
        return trapezoid_from_flat_area(flat_area, flat_time; delay, rise_time,
            fall_time, sys, max_grad, max_slew)
    else
        if flat_time !== nothing
            return trapezoid_from_amplitude(amplitude, flat_time; delay, rise_time,
                fall_time, sys, max_grad, max_slew)
        elseif duration !== nothing
            return trapezoid_from_amplitude_duration(amplitude, duration; delay,
                rise_time, fall_time, sys, max_grad, max_slew)
        else
            error("Supply flat_time or duration with amplitude.")
        end
    end
end

function trapezoid_from_flat_area(flat_area, flat_time; delay=0.0, rise_time=nothing,
    fall_time=nothing, sys=Scanner(), max_grad=sys.Gmax, max_slew=sys.Smax)
    return trapezoid_from_amplitude(
        flat_area / flat_time, flat_time; delay, rise_time, fall_time, sys, max_grad,
        max_slew,
    )
end

function trapezoid_from_amplitude(amplitude, flat_time; delay=0.0, rise_time=nothing,
    fall_time=nothing, sys=Scanner(), max_grad=sys.Gmax, max_slew=sys.Smax)
    if rise_time === nothing
        fall_time !== nothing && error("Supply rise_time when fall_time is specified.")
        rise_time = slew_limited_rise_time(amplitude; sys, max_slew)
    end
    fall_time = something(fall_time, rise_time)
    flat_time < 0 && error("Gradient ramps exceed duration.")
    grad = Grad(amplitude, flat_time, rise_time, fall_time, delay)
    check_hw_limits(grad; max_grad, max_slew=Inf)
    return grad
end

function trapezoid_from_amplitude_duration(amplitude, duration; delay=0.0,
    rise_time=nothing, fall_time=nothing, sys=Scanner(), max_grad=sys.Gmax,
    max_slew=sys.Smax)
    if rise_time === nothing
        fall_time !== nothing && error("Supply rise_time when fall_time is specified.")
        rise_time = slew_limited_rise_time(amplitude; sys, max_slew)
    end
    fall_time = something(fall_time, rise_time)
    flat_time = duration - rise_time - fall_time
    return trapezoid_from_amplitude(
        amplitude, flat_time; delay, rise_time, fall_time, sys, max_grad, max_slew,
    )
end

function trapezoid_from_area(target_area, duration; delay=0.0, sys=Scanner(),
    max_grad=sys.Gmax, max_slew=sys.Smax, rise_time=nothing, fall_time=nothing)
    if rise_time === nothing
        fall_time !== nothing && error("Supply rise_time when fall_time is specified.")
        amplitude = solve_max_slew_amplitude(target_area, duration, max_slew)
        rise_time = slew_limited_rise_time(amplitude; sys, max_slew)
    end
    fall_time = something(fall_time, rise_time)
    flat_time = duration - rise_time - fall_time
    flat_time < 0 && error("Gradient ramps exceed duration.")
    grad = Grad(1.0, flat_time, rise_time, fall_time, delay)
    amplitude = target_area / area(grad)
    grad = amplitude * grad
    check_hw_limits(grad; max_grad, max_slew=Inf)
    return grad
end

function shortest_trapezoid(target_area, delay, raster, max_grad, max_slew)
    rise_time = ramp_time_to_raster(sqrt(abs(target_area) / max_slew), raster)
    amplitude = target_area / rise_time
    effective_time = rise_time
    if abs(amplitude) > max_grad
        effective_time = ceil_to_raster(abs(target_area) / max_grad, raster)
        amplitude = target_area / effective_time
        rise_time = ramp_time_to_raster(abs(amplitude) / max_slew, raster)
    end
    grad = Grad(amplitude, effective_time - rise_time, rise_time, rise_time, delay)
    check_hw_limits(grad; max_grad, max_slew=Inf)
    return grad
end

function solve_max_slew_amplitude(target_area, duration, max_slew)
    c = 1 / max_slew
    disc = duration^2 - 4 * abs(target_area) * c
    disc <= 0 && error("Requested area is too large for this duration.")
    return sign(target_area) * (duration - sqrt(disc)) / (2c)
end

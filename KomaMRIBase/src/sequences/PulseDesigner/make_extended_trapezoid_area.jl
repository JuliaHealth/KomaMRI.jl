"""
    seq = build_extended_trapezoid_area(axis::Symbol, grad_start, grad_end, area; kwargs...)

Return a one-block `Sequence` with the shortest Pulseq-style extended trapezoid
matching the requested edge amplitudes and area. See `make_extended_trapezoid_area`
for gradient design keywords.

# Arguments
- `axis::Symbol`: Gradient axis, one of `:x`, `:y`, or `:z`.
- `grad_start`: Starting gradient amplitude. Plain numbers use Pulseq units. [`Hz/m`]
- `grad_end`: Ending gradient amplitude. Plain numbers use Pulseq units. [`Hz/m`]
- `area`: Target gradient area. Plain numbers use Pulseq units. [`1/m`]

# Returns
- `seq`: Sequence containing the gradient block.
"""
function build_extended_trapezoid_area(
    axis::Symbol, grad_start, grad_end, area; sys=Scanner(),
)
    seq = Sequence(sys)
    axis in (:x, :y, :z) || error("Gradient axis must be :x, :y, or :z.")
    grad = make_extended_trapezoid_area(grad_start, grad_end, area; sys)
    addblock!(seq; (; axis => grad)...)
    return seq
end

"""
    grad = make_extended_trapezoid_area(grad_start, grad_end, area; sys=Scanner())

Return the shortest Pulseq-style extended trapezoid matching the requested
edge amplitudes and area.

# Arguments
- `grad_start`: Starting gradient amplitude. Plain numbers use Pulseq units. [`Hz/m`]
- `grad_end`: Ending gradient amplitude. Plain numbers use Pulseq units. [`Hz/m`]
- `area`: Target gradient area. Plain numbers use Pulseq units. [`1/m`]

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `grad`: Extended trapezoid gradient event.
"""
function make_extended_trapezoid_area(grad_start, grad_end, target_area; sys=Scanner())
    grad_start  = to_SI(grad_start, PulseqUnitsDefault())
    grad_end    = to_SI(grad_end, PulseqUnitsDefault())
    target_area = to_SI(target_area, PulseqUnitsDefault())
    return extended_trapezoid_area(grad_start, grad_end, target_area; sys)
end

function extended_trapezoid_area(grad_start, grad_end, target_area; sys=Scanner())
    if iszero(grad_start) && iszero(grad_end) && iszero(target_area)
        return Grad(0.0, 0.0)
    end
    # MATLAB Pulseq uses a 1% margin here when searching for a feasible solution.
    max_slew = 0.99 * sys.Smax
    max_grad = 0.99 * sys.Gmax
    raster = sys.GR_Δt
    min_duration = max(
        slew_limited_ramp_samples(grad_end - grad_start, max_slew, raster), 2,
    )
    max_duration = maximum((
        slew_limited_ramp_samples(grad_start, max_slew, raster),
        slew_limited_ramp_samples(grad_end, max_slew, raster),
        min_duration,
    ))
    solution = nothing
    for duration in min_duration:max_duration
        solution = extended_trapezoid_area_solution(
            duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
        )
        solution === nothing || break
    end
    if solution === nothing
        duration = max_duration
        while solution === nothing
            duration *= 2
            solution = extended_trapezoid_area_solution(
                duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
            )
        end
        solution = extended_trapezoid_area_binary_search(
            duration ÷ 2, duration, target_area, grad_start, grad_end,
            max_slew, max_grad, raster,
        )
    end
    ramp_up, flat, ramp_down, amplitude = solution
    times, amplitudes = extended_trapezoid_area_points(
        ramp_up, flat, ramp_down, grad_start, amplitude, grad_end, raster,
    )
    grad = extended_trapezoid(times, amplitudes; sys)
    # Only allow floating-point roundoff here; Pulseq's coarse 1e-3 check is too loose.
    isapprox(area(grad), target_area; rtol=sqrt(eps(Float64)), atol=eps(Float64)) ||
        error("Could not find a solution.")
    return grad
end

function extended_trapezoid_area_binary_search(
    low, high, target_area, grad_start, grad_end, max_slew, max_grad, raster,
)
    while low < high - 1
        mid = (low + high) ÷ 2
        solution = extended_trapezoid_area_solution(
            mid, target_area, grad_start, grad_end, max_slew, max_grad, raster,
        )
        solution === nothing ? low = mid : high = mid
    end
    return extended_trapezoid_area_solution(
        high, target_area, grad_start, grad_end, max_slew, max_grad, raster,
    )
end

function extended_trapezoid_area_points(
    ramp_up, flat, ramp_down, grad_start, amplitude, grad_end, raster,
)
    if flat > 0
        times = cumsum([0.0, ramp_up, flat, ramp_down] .* raster)
        amplitudes = [grad_start, amplitude, amplitude, grad_end]
    else
        times = cumsum([0.0, ramp_up, ramp_down] .* raster)
        amplitudes = [grad_start, amplitude, grad_end]
    end
    return times, amplitudes
end

function extended_trapezoid_area_solution(
    duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
)
    extended_trapezoid_area_is_reachable(
        duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
    ) || return nothing
    solution = extended_trapezoid_area_max_grad_solution(
        duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
    )
    solution === nothing || return solution
    return extended_trapezoid_area_best_solution(
        duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
    )
end

function extended_trapezoid_area_is_reachable(
    duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
)
    iszero(target_area) && return true
    area_sign = sign(target_area)
    amplitude = area_sign * max_grad
    ramp_up = abs(amplitude - grad_start) / max_slew / raster
    ramp_down = abs(amplitude - grad_end) / max_slew / raster
    flat = max(duration - ramp_up - ramp_down, 0)
    max_area = ramp_up * (amplitude + grad_start) +
        ramp_down * (amplitude + grad_end) + 2 * flat * amplitude
    return abs(2 * target_area / raster) <= abs(max_area)
end

function extended_trapezoid_area_max_grad_solution(
    duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
)
    iszero(target_area) && return nothing
    area_sign = sign(target_area)
    ramp_up = (duration * max_slew * raster + area_sign * (grad_end - grad_start)) /
        (2 * max_slew * raster)
    peak = area_sign * grad_start + ramp_up * max_slew * raster
    (peak <= max_grad || isapprox(peak, max_grad; rtol=sqrt(eps(Float64)), atol=0)) &&
        return nothing
    ramp_up = round(Int, abs(grad_start - area_sign * max_grad) / max_slew / raster)
    ramp_down = round(Int, abs(grad_end - area_sign * max_grad) / max_slew / raster)
    flat = duration - ramp_up - ramp_down
    flat > 0 || return nothing
    amplitude = -(ramp_up * raster * grad_start + ramp_down * raster * grad_end -
        2 * target_area) / ((ramp_up + 2 * flat + ramp_down) * raster)
    times, amplitudes = extended_trapezoid_area_points(
        ramp_up, flat, ramp_down, grad_start, amplitude, grad_end, raster,
    )
    slew = diff(amplitudes) ./ diff(times)
    slew_peak = maximum(abs, slew)
    amplitude_peak = maximum(abs, amplitudes)
    (slew_peak <= max_slew || isapprox(slew_peak, max_slew; rtol=sqrt(eps(Float64)), atol=0)) &&
        (amplitude_peak <= max_grad ||
            isapprox(amplitude_peak, max_grad; rtol=sqrt(eps(Float64)), atol=0)) ||
        return nothing
    return ramp_up, flat, ramp_down, amplitude
end

function extended_trapezoid_area_best_solution(
    duration, target_area, grad_start, grad_end, max_slew, max_grad, raster,
)
    ramp_up_limit =
        ceil(Int, max(abs(max_grad - grad_start), abs(-max_grad - grad_start)) /
            max_slew / raster) + 1
    ramp_down_limit =
        ceil(Int, max(abs(max_grad - grad_end), abs(-max_grad - grad_end)) /
            max_slew / raster) + 1
    best = nothing
    best_slew = Inf
    for ramp_up in 1:ramp_up_limit, ramp_down in 1:ramp_down_limit
        flat = duration - ramp_up - ramp_down
        flat > 0 && ramp_up > 0 && ramp_down > 0 || continue
        candidate = extended_trapezoid_area_candidate(
            ramp_up, flat, ramp_down, target_area, grad_start, grad_end, raster,
            max_slew, max_grad,
        )
        candidate === nothing && continue
        slew, amplitude = candidate
        slew < best_slew || continue
        best = ramp_up, flat, ramp_down, amplitude
        best_slew = slew
    end
    for ramp_up in 1:(duration - 1)
        ramp_down = duration - ramp_up
        candidate = extended_trapezoid_area_candidate(
            ramp_up, 0, ramp_down, target_area, grad_start, grad_end, raster,
            max_slew, max_grad,
        )
        candidate === nothing && continue
        slew, amplitude = candidate
        best !== nothing &&
            isapprox(slew, best_slew; rtol=sqrt(eps(Float64)), atol=0) && continue
        slew < best_slew || continue
        best = ramp_up, 0, ramp_down, amplitude
        best_slew = slew
    end
    return best
end

function extended_trapezoid_area_candidate(
    ramp_up, flat, ramp_down, target_area, grad_start, grad_end, raster,
    max_slew, max_grad,
)
    amplitude = -(ramp_up * raster * grad_start + ramp_down * raster * grad_end -
        2 * target_area) / ((ramp_up + 2 * flat + ramp_down) * raster)
    slew_up = abs(grad_start - amplitude) / (ramp_up * raster)
    slew_down = abs(grad_end - amplitude) / (ramp_down * raster)
    amplitude_peak = abs(amplitude)
    rtol = sqrt(eps(Float64))
    (amplitude_peak <= max_grad || isapprox(amplitude_peak, max_grad; rtol, atol=0)) ||
        return nothing
    (slew_up <= max_slew || isapprox(slew_up, max_slew; rtol, atol=0)) &&
        (slew_down <= max_slew || isapprox(slew_down, max_slew; rtol, atol=0)) ||
        return nothing
    return slew_up + slew_down, amplitude
end

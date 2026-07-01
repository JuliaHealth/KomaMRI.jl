"""
    check_hw_limits(seq)
    check_hw_limits(seq, sys)

Checks event-local hardware limits:
- maximum RF amplitude `|B1|`
- maximum gradient amplitude `|G|`
- maximum gradient slew rate `|dG/dt|`

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `sys`: (`::Scanner`) Scanner struct. If omitted, hardware metadata stored in
  `seq.DEF` is used.
"""
function check_hw_limits(seq::Sequence)
    return check_hw_limits(seq, _sequence_scanner_from_def(seq.DEF))
end

function check_hw_limits(seq::Sequence, sys::Scanner)
    return check_hw_limits(seq.GR, seq.RF, seq.ADC, sys)
end

_absmax(x::Number) = abs(x)
_absmax(x::AbstractArray{<:Number}) = maximum(abs, x)

function check_hw_limits(gr::Grad, sys::Scanner)
    return check_hw_limits(gr; max_grad=sys.Gmax, max_slew=sys.Smax)
end

function check_hw_limits(gr::Grad; max_grad=Inf, max_slew=Inf, name="Gradient")
    is_GR_on(gr) || return nothing
    rtol = sqrt(eps(Float64))
    if isfinite(max_grad)
        grad_peak = max(abs(gr.first), _absmax(gr.A), abs(gr.last))
        if grad_peak > max_grad && !isapprox(grad_peak, max_grad; rtol, atol=0)
            limit = max_grad * 1e3
            error("$name amplitude is greater than the maximum gradient amplitude of the scanner ($limit mT/m).")
        end
    end
    if isfinite(max_slew)
        grad_slew = _max_gradient_slew(gr)
        if grad_slew > max_slew && !isapprox(grad_slew, max_slew; rtol, atol=0)
            error("$name slew rate is greater than the maximum gradient slew rate of the scanner ($(max_slew) mT/m/ms).")
        end
    end
    return nothing
end

function check_hw_limits(gr::AbstractMatrix{G}, rf::AbstractMatrix{R}, ::AbstractVector{A}, sys::Scanner) where {G<:Grad,R<:RF,A<:ADC}
    check_rf = isfinite(sys.B1)
    check_grad = isfinite(sys.Gmax)
    check_slew = isfinite(sys.Smax)
    (check_rf || check_grad || check_slew) || return nothing
    axis_names = ("x", "y", "z")
    rtol = sqrt(eps(Float64))
    for i in axes(rf, 2)
        if check_rf
            rf_i = rf[1, i]
            rf_peak = is_RF_on(rf_i) ? _absmax(rf_i.A) : 0.0
            if rf_peak > sys.B1 && !isapprox(rf_peak, sys.B1; rtol, atol=0)
                error("RF amplitude for block $i is greater than the maximum RF amplitude of the scanner ($(sys.B1 * 1e6) μT).")
            end
        end
        if check_grad || check_slew
            for gi in 1:3
                gr_i = gr[gi, i]
                is_GR_on(gr_i) || continue
                max_grad = check_grad ? sys.Gmax : Inf
                max_slew = check_slew ? sys.Smax : Inf
                name = "$(axis_names[gi]) gradient for block $i"
                check_hw_limits(gr_i; max_grad, max_slew, name)
            end
        end
    end
    return nothing
end

function _max_gradient_slew(gr::TrapezoidalGrad)
    max_slew = gr.rise > 0 ? abs(gr.A - gr.first) / gr.rise : 0.0
    if gr.fall > 0
        max_slew = max(max_slew, abs(gr.last - gr.A) / gr.fall)
    end
    return max_slew
end

function _max_gradient_slew(gr::UniformlySampledGrad)
    n = length(gr.A)
    n == 0 && return 0.0
    thresh = max(100 * MIN_RISE_TIME, sqrt(eps(Float64)))
    max_slew = gr.rise > 0 ? abs(gr.A[1] - gr.first) / gr.rise : 0.0
    if n > 1
        Δt = gr.T / (n - 1)
        @inbounds for i in 1:(n - 1)
            abs(Δt) > thresh || continue
            max_slew = max(max_slew, abs((gr.A[i + 1] - gr.A[i]) / Δt))
        end
    end
    if gr.fall > 0
        max_slew = max(max_slew, abs(gr.last - gr.A[end]) / gr.fall)
    end
    return max_slew
end

function _max_gradient_slew(gr::TimeShapedGrad)
    n = length(gr.A)
    n == 0 && return 0.0
    thresh = max(100 * MIN_RISE_TIME, sqrt(eps(Float64)))
    max_slew = gr.rise > 0 ? abs(gr.A[1] - gr.first) / gr.rise : 0.0
    @inbounds for i in 1:min(length(gr.T), n - 1)
        abs(gr.T[i]) > thresh || continue
        max_slew = max(max_slew, abs((gr.A[i + 1] - gr.A[i]) / gr.T[i]))
    end
    if gr.fall > 0
        max_slew = max(max_slew, abs(gr.last - gr.A[end]) / gr.fall)
    end
    return max_slew
end

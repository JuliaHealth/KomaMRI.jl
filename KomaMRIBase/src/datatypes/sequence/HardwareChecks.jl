"""
    check_hw_limits(seq)
    check_hw_limits(seq, sys)

Checks event-local hardware limits:
- maximum RF amplitude `|B1|`
- maximum gradient amplitude `|G|`
- maximum gradient slew rate `|dG/dt|`
- minimum ADC dwell time

It does not enforce RF/ADC dead times or RF ring-down.

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

function check_hw_limits(gr::AbstractMatrix{G}, rf::AbstractMatrix{R}, adc::AbstractVector{ADC}, sys::Scanner) where {G<:Grad,R<:RF}
    check_rf = isfinite(sys.B1)
    check_grad = isfinite(sys.Gmax)
    check_slew = isfinite(sys.Smax)
    check_adc = sys.ADC_Δt > 0
    (check_rf || check_grad || check_slew || check_adc) || return nothing
    axis_names = ("x", "y", "z")
    rtol = sqrt(eps(Float64))
    for i in eachindex(adc)
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
                if check_grad
                    grad_peak = max(abs(gr_i.first), _absmax(gr_i.A), abs(gr_i.last))
                    (grad_peak <= sys.Gmax || isapprox(grad_peak, sys.Gmax; rtol, atol=0)) || error("$(axis_names[gi]) gradient amplitude for block $i is greater than the maximum gradient amplitude of the scanner ($(sys.Gmax * 1e3) mT/m).")
                end
                if check_slew
                    grad_slew = _max_gradient_slew(gr_i)
                    (grad_slew <= sys.Smax || isapprox(grad_slew, sys.Smax; rtol, atol=0)) || error("$(axis_names[gi]) gradient slew rate for block $i is greater than the maximum gradient slew rate of the scanner ($(sys.Smax) mT/m/ms).")
                end
            end
        end
        if check_adc
            adc_i = adc[i]
            is_ADC_on(adc_i) || continue
            dwell = adc_i.N == 1 ? adc_i.T : adc_i.T / (adc_i.N - 1)
            if dwell < sys.ADC_Δt
                error("ADC dwell time $(dwell * 1e6) μs for block $i is less than the minimum ADC dwell time of the scanner ($(sys.ADC_Δt * 1e6) μs).")
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

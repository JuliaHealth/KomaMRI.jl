"""
    check_hw_limits(seq)
    check_hw_limits(seq, sys)

Checks event-local hardware limits:
- maximum RF amplitude `|B1|`
- maximum gradient amplitude `|G|`
- maximum gradient slew rate `|dG/dt|`
- ADC sample count accepted by scanners

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `sys`: (`::Scanner`) Scanner struct. If omitted, hardware metadata stored in
  `seq.DEF` is used.
"""
function check_hw_limits(seq::Sequence)
    return check_hw_limits(seq, _sequence_scanner_from_def(seq.DEF))
end

function check_hw_limits(seq::Sequence, sys::Scanner)
    check_rf = isfinite(sys.B1)
    check_grad = isfinite(sys.Gmax)
    check_slew = isfinite(sys.Smax)
    check_fields = check_rf || check_grad || check_slew
    axis_names = ("x", "y", "z")
    for i in 1:length(seq)
        adc_i = seq.ADC[i]
        if is_ADC_on(adc_i)
            isodd(adc_i.N) && error("ADC for block $i has an odd number of samples ($(adc_i.N)).")
        end
        if check_fields
            if check_rf
                rf_i = seq.RF[1, i]
                rf_peak = is_RF_on(rf_i) ? maximum(abs, ampls(rf_i)) : 0.0
                if rf_peak > sys.B1 && !(rf_peak ≈ sys.B1)
                    error("RF amplitude for block $i is greater than the maximum RF amplitude of the scanner ($(sys.B1 * 1e6) μT).")
                end
            end
            if check_grad || check_slew
                for gi in 1:size(seq.GR, 1)
                    gr_i = seq.GR[gi, i]
                    is_GR_on(gr_i) || continue
                    max_grad = check_grad ? sys.Gmax : Inf
                    max_slew = check_slew ? sys.Smax : Inf
                    name = "$(axis_names[gi]) gradient for block $i"
                    check_hw_limits(gr_i; max_grad, max_slew, name)
                end
            end
        end
    end
    return nothing
end

function check_hw_limits(gr::Grad; max_grad=Inf, max_slew=Inf, name="Gradient")
    is_GR_on(gr) || return nothing
    if isfinite(max_grad)
        grad_peak = maximum(abs, ampls(gr))
        if grad_peak > max_grad && !(grad_peak ≈ max_grad)
            limit = max_grad * 1e3
            error("$name amplitude is greater than the maximum gradient amplitude of the scanner ($limit mT/m).")
        end
    end
    if isfinite(max_slew)
        grad_slew = _max_gradient_slew(gr)
        if grad_slew > max_slew && !(grad_slew ≈ max_slew)
            error("$name slew rate is greater than the maximum gradient slew rate of the scanner ($(max_slew) mT/m/ms).")
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

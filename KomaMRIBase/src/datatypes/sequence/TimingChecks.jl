const PULSEQ_TIME_TOL = 1e-12
const PULSEQ_DIVISION_TOL = 1e-9
const PULSEQ_ADC_DWELL_TOL = 1e-10

function _sequence_raster_from_def(def)
    missing = [key for key in PULSEQ_RASTER_DEFINITION_KEYS if !haskey(def, key)]
    isempty(missing) || error("Sequence has no Pulseq raster definitions. Use `Sequence()` or `Sequence(sys)`.")
    return (;
        BlockDurationRaster=def["BlockDurationRaster"],
        GradientRasterTime=def["GradientRasterTime"],
        RadiofrequencyRasterTime=def["RadiofrequencyRasterTime"],
        AdcRasterTime=def["AdcRasterTime"],
    )
end

_sequence_raster_from_sys(sys::Scanner) = (;
    BlockDurationRaster=sys.DUR_Δt,
    GradientRasterTime=sys.GR_Δt,
    RadiofrequencyRasterTime=sys.RF_Δt,
    AdcRasterTime=sys.ADC_Δt,
)

function _check_raster_multiple(t, raster, block_id, label; tol=PULSEQ_DIVISION_TOL)
    raster > 0 || error("Raster time must be positive.")
    abs(t / raster - round(t / raster)) <= tol && return nothing
    error("Block $block_id $label ($(t) s) is not aligned to raster $(raster) s.")
end

_compact_grad_timing(start, step, n, raster) =
    isapprox(start, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) &&
    (isapprox(step, raster; rtol=0, atol=PULSEQ_TIME_TOL) ||
     (isodd(n) && isapprox(step, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL)))

function _check_uniform_timing(start, step, n, raster, block_id, label)
    n == 0 && return nothing
    _check_raster_multiple(start, raster, block_id, label)
    n == 1 && return nothing
    _check_raster_multiple(step, raster, block_id, label)
    return nothing
end

function _check_cumulative_timing(T, n, start, raster, block_id, label)
    n == 0 && return nothing
    lenT = length(T)
    (lenT == n - 1 || lenT == n) || throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $lenT."))
    t = start
    _check_raster_multiple(t, raster, block_id, label)
    for i in 1:(n - 1)
        t += T[i]
        _check_raster_multiple(t, raster, block_id, label)
    end
    return nothing
end

function _is_cumulative_compact_grad(T, n, start, raster)
    n == 0 && return false
    lenT = length(T)
    (lenT == n - 1 || lenT == n) || throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $lenT."))
    step = n > 1 ? T[1] : zero(start)
    _compact_grad_timing(start, step, n, raster) || return false
    t = start
    for i in 1:(n - 1)
        t += T[i]
        target = isapprox(step, raster; rtol=0, atol=PULSEQ_TIME_TOL) ? (i + 0.5) * raster : (i + 1) * (raster / 2)
        isapprox(t, target; rtol=0, atol=PULSEQ_TIME_TOL) || return false
    end
    return true
end

function _check_grad_timing(gr::TrapezoidalGrad, raster, block_id, name)
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, raster, block_id, "$name-gradient rise time")
    _check_raster_multiple(gr.T, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, raster, block_id, "$name-gradient fall time")
    return nothing
end

function _check_grad_timing(gr::UniformlySampledGrad, raster, block_id, name)
    n = length(gr.A)
    step = n > 1 ? gr.T / (n - 1) : zero(gr.T)
    compact = _compact_grad_timing(gr.rise, step, n, raster) &&
        isapprox(gr.fall, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL)
    edge_raster = compact ? raster / 2 : raster
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, edge_raster, block_id, "$name-gradient rise time")
    compact || _check_uniform_timing(gr.rise, step, n, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, edge_raster, block_id, "$name-gradient fall time")
    return nothing
end

function _check_grad_timing(gr::TimeShapedGrad, raster, block_id, name)
    n = length(gr.A)
    compact = isapprox(gr.fall, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) &&
        _is_cumulative_compact_grad(gr.T, n, gr.rise, raster)
    edge_raster = compact ? raster / 2 : raster
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, edge_raster, block_id, "$name-gradient rise time")
    compact || _check_cumulative_timing(gr.T, n, gr.rise, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, edge_raster, block_id, "$name-gradient fall time")
    return nothing
end

function _check_rf_timing(rf::BlockPulseRF, raster, block_id)
    _check_raster_multiple(rf.delay, raster, block_id, "RF delay")
    _check_raster_multiple(rf.T, raster, block_id, "RF timing")
    return nothing
end

function _check_rf_timing(rf::UniformlySampledRF, raster, block_id)
    n = length(rf.A)
    step = n > 1 ? rf.T / (n - 1) : zero(rf.T)
    compact = isapprox(step, raster; rtol=0, atol=PULSEQ_TIME_TOL) ||
        (isodd(n) && isapprox(step, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL))
    delay = compact ? rf.delay - raster / 2 : rf.delay
    _check_raster_multiple(delay, raster, block_id, "RF delay")
    compact || _check_uniform_timing(0.0, step, n, raster, block_id, "RF timing")
    return nothing
end

function _check_rf_timing(rf::TimeShapedRF, raster, block_id)
    _check_raster_multiple(rf.delay, raster, block_id, "RF delay")
    _check_cumulative_timing(rf.T, length(rf.A), 0.0, raster, block_id, "RF timing")
    return nothing
end

function _check_adc_dwell(dwell, raster, block_id)
    raster > 0 || error("Raster time must be positive.")
    dwell + PULSEQ_TIME_TOL >= raster || error("Block $block_id ADC dwell ($(dwell) s) is shorter than raster $(raster) s.")
    rounded = round(dwell / raster) * raster
    abs(rounded - dwell) <= PULSEQ_ADC_DWELL_TOL && return nothing
    error("Block $block_id ADC dwell ($(dwell) s) is not aligned to raster $(raster) s.")
end

function _check_adc_timing(adc, adc_raster, rf_raster, block_id)
    dwell = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
    _check_adc_dwell(dwell, adc_raster, block_id)
    adc.delay + PULSEQ_TIME_TOL >= dwell / 2 || error("Block $block_id ADC delay ($(adc.delay) s) is smaller than dwell/2 ($(dwell / 2) s).")
    _check_raster_multiple(adc.delay - dwell / 2, rf_raster, block_id, "ADC delay")
    return nothing
end

"""
    check_timing(seq)
    check_timing(seq, sys)

Check Pulseq timing alignment. The one-argument method uses raster definitions
stored in `seq.DEF`; the two-argument method uses raster definitions from `sys`.
"""
function check_timing(seq::Sequence)
    return check_timing(seq, _sequence_raster_from_def(seq.DEF))
end

function check_timing(seq::Sequence, sys::Scanner)
    return check_timing(seq, _sequence_raster_from_sys(sys))
end

function check_timing(seq::Sequence, raster::NamedTuple)
    for i in eachindex(seq.DUR)
        _check_raster_multiple(seq.DUR[i], raster.BlockDurationRaster, i, "duration")
        rf = seq.RF[1, i]
        if is_RF_on(rf)
            _check_rf_timing(rf, raster.RadiofrequencyRasterTime, i)
        end
        axis_names = ("x", "y", "z")
        for axis in axes(seq.GR, 1)
            gr = seq.GR[axis, i]
            is_GR_on(gr) || continue
            name = axis <= length(axis_names) ? axis_names[axis] : string(axis)
            _check_grad_timing(gr, raster.GradientRasterTime, i, name)
        end
        adc = seq.ADC[i]
        is_ADC_on(adc) && _check_adc_timing(adc, raster.AdcRasterTime, raster.RadiofrequencyRasterTime, i)
    end
    return nothing
end

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

_is_raster_multiple(t, raster; tol=PULSEQ_DIVISION_TOL) = abs(t / raster - round(t / raster)) <= tol

function _check_raster_multiple(t, raster, message; tol=PULSEQ_DIVISION_TOL)
    raster > 0 || error("Raster time must be positive.")
    _is_raster_multiple(t, raster; tol) && return nothing
    error("$message ($(t) s) is not aligned to raster $(raster) s.")
end

function _check_raster_multiple(ts::AbstractVector, raster, message; tol=PULSEQ_DIVISION_TOL)
    for t in ts
        _check_raster_multiple(t, raster, message; tol)
    end
    return nothing
end

_grad_sample_times(gr) = _shape_times(gr.A, gr.T)

function _grad_sample_times(gr::UniformlySampledGrad)
    n = length(gr.A)
    interval = n > 1 ? gr.T / (n - 1) : 0.0
    return gr.rise .+ (0:(n - 1)) .* interval
end

_grad_sample_times(gr::TimeShapedGrad) = gr.rise .+ cumsum(vcat(0.0, gr.T))

function _grad_is_compact_for_timing(gr, tt, raster)
    isapprox(gr.rise, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return false
    isapprox(gr.fall, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return false
    n = length(tt)
    default_tt = ((1:n) .- 0.5) .* raster
    oversampled_tt = (1:n) .* (raster / 2)
    return isapprox(tt, default_tt; rtol=0, atol=PULSEQ_TIME_TOL) ||
        (isodd(n) && isapprox(tt, oversampled_tt; rtol=0, atol=PULSEQ_TIME_TOL))
end

function _check_grad_timing(gr, raster, block_id, name)
    tt = _grad_sample_times(gr)
    compact = _grad_is_compact_for_timing(gr, tt, raster)
    edge_raster = compact ? raster / 2 : raster
    _check_raster_multiple(gr.delay, raster, "Block $block_id $name-gradient delay")
    _check_raster_multiple(gr.rise, edge_raster, "Block $block_id $name-gradient rise time")
    compact || _check_raster_multiple(tt, raster, "Block $block_id $name-gradient timing")
    _check_raster_multiple(gr.fall, edge_raster, "Block $block_id $name-gradient fall time")
    return nothing
end

_rf_is_compact_for_timing(::TimeShapedRF, tt, raster) = false
function _rf_is_compact_for_timing(rf, tt, raster)
    default_tt = (0:(length(tt) - 1)) .* raster
    oversampled_tt = (0:(length(tt) - 1)) .* (raster / 2)
    return isapprox(tt, default_tt; rtol=0, atol=PULSEQ_TIME_TOL) ||
        (isodd(length(tt)) && isapprox(tt, oversampled_tt; rtol=0, atol=PULSEQ_TIME_TOL))
end

function _check_rf_timing(rf, raster, block_id)
    tt = _shape_times(rf.A, rf.T)
    compact = _rf_is_compact_for_timing(rf, tt, raster)
    delay = compact ? rf.delay - raster / 2 : rf.delay
    _check_raster_multiple(delay, raster, "Block $block_id RF delay")
    compact || _check_raster_multiple(tt, raster, "Block $block_id RF timing")
    return nothing
end

function _check_adc_dwell(dwell, raster, block_id)
    raster > 0 || error("Raster time must be positive.")
    dwell + PULSEQ_TIME_TOL >= raster || error("Block $block_id ADC dwell ($(dwell) s) is shorter than raster $(raster) s.")
    rounded = round(dwell / raster) * raster
    abs(rounded - dwell) <= PULSEQ_ADC_DWELL_TOL && return nothing
    error("Block $block_id ADC dwell ($(dwell) s) is not aligned to raster $(raster) s.")
end

function _check_adc_timing(adc, raster, rf_raster, block_id)
    dwell = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
    _check_adc_dwell(dwell, raster, block_id)
    adc.delay + PULSEQ_TIME_TOL >= dwell / 2 || error("Block $block_id ADC delay ($(adc.delay) s) is smaller than dwell/2 ($(dwell / 2) s).")
    _check_raster_multiple(adc.delay - dwell / 2, rf_raster, "Block $block_id ADC delay")
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
        _check_raster_multiple(seq.DUR[i], raster.BlockDurationRaster, "Block $i duration")
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

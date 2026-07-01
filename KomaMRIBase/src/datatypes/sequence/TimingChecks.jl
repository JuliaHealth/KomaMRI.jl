_sequence_timing_from_sys(sys::Scanner) = (;
    BlockDurationRaster=sys.DUR_Δt,
    GradientRasterTime=sys.GR_Δt,
    RadiofrequencyRasterTime=sys.RF_Δt,
    AdcRasterTime=sys.ADC_Δt,
    RfRingdownTime=sys.RF_ring_down_time,
    RfDeadTime=sys.RF_dead_time,
    AdcDeadTime=sys.ADC_dead_time,
)

function _check_block_fit(event_end, block_duration, block_id, label)
    event_end <= block_duration + PULSEQ_TIME_TOL && return nothing
    error("Block $block_id $label exceeds the block duration.")
end

function _check_raster_multiple(t, raster, block_id, label; tol=PULSEQ_DIVISION_TOL)
    raster > 0 || error("Raster time must be positive.")
    abs(t / raster - round(t / raster)) <= tol && return nothing
    error("Block $block_id $label ($(t) s) is not aligned to raster $(raster) s.")
end

function _check_sample_timing(start::Number, step::Number, n, raster, block_id, label)
    n == 0 && return nothing
    _check_raster_multiple(start, raster, block_id, label)
    n == 1 && return nothing
    _check_raster_multiple(step, raster, block_id, label)
    return nothing
end

function _check_sample_timing(T::AbstractVector, n, start, raster, block_id, label)
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

function _check_timing(gr::TrapezoidalGrad, timing, block_id, name, block_duration)
    raster = timing.GradientRasterTime
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, raster, block_id, "$name-gradient rise time")
    _check_raster_multiple(gr.T, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, raster, block_id, "$name-gradient fall time")
    _check_block_fit(dur(gr), block_duration, block_id, "$name-gradient event")
    return nothing
end

function _check_timing(gr::UniformlySampledGrad, timing, block_id, name, block_duration)
    raster = timing.GradientRasterTime
    n = length(gr.A)
    step = n > 1 ? gr.T / (n - 1) : zero(gr.T)
    compact_interval = compact_sample_interval(gr, raster)
    edge_raster = isnothing(compact_interval) ? raster : raster / 2
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, edge_raster, block_id, "$name-gradient rise time")
    isnothing(compact_interval) && _check_sample_timing(gr.rise, step, n, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, edge_raster, block_id, "$name-gradient fall time")
    _check_block_fit(dur(gr), block_duration, block_id, "$name-gradient event")
    return nothing
end

function _check_timing(gr::TimeShapedGrad, timing, block_id, name, block_duration)
    raster = timing.GradientRasterTime
    n = length(gr.A)
    compact_interval = compact_sample_interval(gr, raster)
    edge_raster = isnothing(compact_interval) ? raster : raster / 2
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, edge_raster, block_id, "$name-gradient rise time")
    isnothing(compact_interval) && _check_sample_timing(gr.T, n, gr.rise, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, edge_raster, block_id, "$name-gradient fall time")
    _check_block_fit(dur(gr), block_duration, block_id, "$name-gradient event")
    return nothing
end

function _check_timing(rf::BlockPulseRF, timing, sys, block_id, block_duration)
    raster = timing.RadiofrequencyRasterTime
    _check_raster_multiple(delay(rf, sys), raster, block_id, "RF delay")
    _check_raster_multiple(sum(rf.T), raster, block_id, "RF timing")
    _check_rf_duration(rf, sys, block_id, block_duration)
    return nothing
end

function _check_timing(rf::UniformlySampledRF, timing, sys, block_id, block_duration)
    raster = timing.RadiofrequencyRasterTime
    n = length(rf.A)
    step = n > 1 ? rf.T / (n - 1) : zero(rf.T)
    compact_interval = compact_sample_interval(rf, raster)
    _check_raster_multiple(delay(rf, sys), raster, block_id, "RF delay")
    isnothing(compact_interval) && _check_sample_timing(0.0, step, n, raster, block_id, "RF timing")
    _check_rf_duration(rf, sys, block_id, block_duration)
    return nothing
end

function _check_timing(rf::TimeShapedRF, timing, sys, block_id, block_duration)
    raster = timing.RadiofrequencyRasterTime
    _check_raster_multiple(delay(rf, sys), raster, block_id, "RF delay")
    _check_rf_duration(rf, sys, block_id, block_duration)
    return nothing
end

function _check_rf_duration(rf, sys, block_id, block_duration)
    _check_block_fit(dur(rf, sys), block_duration, block_id, "RF event")
    if sys.RF_dead_time > 0
        rf_delay = delay(rf, sys)
        rf_delay + PULSEQ_TIME_TOL >= sys.RF_dead_time || error("Block $block_id RF delay ($(rf_delay) s) is smaller than RF dead time $(sys.RF_dead_time) s.")
    end
    return nothing
end

function _check_timing(adc::ADC, timing, sys, block_id, block_duration)
    Δt = dwell(adc, sys)
    adc_raster = timing.AdcRasterTime
    adc_raster > 0 || error("Raster time must be positive.")
    Δt + PULSEQ_TIME_TOL >= adc_raster || error("Block $block_id ADC dwell ($(Δt) s) is shorter than raster $(adc_raster) s.")
    rounded_dwell = round(Δt / adc_raster) * adc_raster
    abs(rounded_dwell - Δt) <= PULSEQ_ADC_DWELL_TOL || error("Block $block_id ADC dwell ($(Δt) s) is not aligned to raster $(adc_raster) s.")
    adc.delay + PULSEQ_TIME_TOL >= Δt / 2 || error("Block $block_id ADC delay ($(adc.delay) s) is smaller than dwell/2 ($(Δt / 2) s).")
    adc_delay = delay(adc, sys)
    _check_raster_multiple(adc_delay, timing.RadiofrequencyRasterTime, block_id, "ADC delay")
    _check_block_fit(dur(adc, sys), block_duration, block_id, "ADC event")
    if sys.ADC_dead_time > 0
        adc_delay + PULSEQ_TIME_TOL >= sys.ADC_dead_time || error("Block $block_id ADC delay ($(adc_delay) s) is smaller than ADC dead time $(sys.ADC_dead_time) s.")
    end
    return nothing
end

function _check_timing(exts::AbstractVector{<:Extension}, block_id, block_duration)
    count(ext -> ext isa Trigger, exts) <= 1 || error("Block $block_id has more than one Trigger extension.")
    count(ext -> ext isa QuaternionRot, exts) <= 1 || error("Block $block_id has more than one ROTATIONS extension.")
    for ext in exts
        _check_block_fit(dur(ext), block_duration, block_id, "$(typeof(ext)) extension")
    end
    return nothing
end

"""
    check_timing(seq)
    check_timing(seq, sys)

Check Pulseq timing alignment, RF/ADC dead times, and RF ring-down. The
one-argument method uses definitions stored in `seq.DEF`; the two-argument
method uses timing definitions from `sys`.
"""
function check_timing(seq::Sequence)
    sys = _sequence_scanner_from_def(seq.DEF)
    return check_timing(seq, _sequence_timing_from_sys(sys), sys)
end

function check_timing(seq::Sequence, sys::Scanner)
    return check_timing(seq, _sequence_timing_from_sys(sys), sys)
end

function check_timing(seq::Sequence, raster::NamedTuple)
    default = Scanner()
    sys = Scanner(
        B0=default.B0,
        B1=default.B1,
        Gmax=default.Gmax,
        Smax=default.Smax,
        ADC_Δt=raster.AdcRasterTime,
        DUR_Δt=raster.BlockDurationRaster,
        GR_Δt=raster.GradientRasterTime,
        RF_Δt=raster.RadiofrequencyRasterTime,
        RF_ring_down_time=get(raster, :RfRingdownTime, 0.0),
        RF_dead_time=get(raster, :RfDeadTime, 0.0),
        ADC_dead_time=get(raster, :AdcDeadTime, 0.0),
    )
    return check_timing(seq, raster, sys)
end

function check_timing(seq::Sequence, raster::NamedTuple, sys::Scanner)
    for i in eachindex(seq.DUR)
        _check_raster_multiple(seq.DUR[i], raster.BlockDurationRaster, i, "duration")
        rf = seq.RF[1, i]
        if is_RF_on(rf)
            _check_timing(rf, raster, sys, i, seq.DUR[i])
        end
        axis_names = ("x", "y", "z")
        for axis in axes(seq.GR, 1)
            gr = seq.GR[axis, i]
            is_GR_on(gr) || continue
            name = axis <= length(axis_names) ? axis_names[axis] : string(axis)
            _check_timing(gr, raster, i, name, seq.DUR[i])
        end
        adc = seq.ADC[i]
        is_ADC_on(adc) && _check_timing(adc, raster, sys, i, seq.DUR[i])
        _check_timing(seq.EXT[i], i, seq.DUR[i])
    end
    return nothing
end

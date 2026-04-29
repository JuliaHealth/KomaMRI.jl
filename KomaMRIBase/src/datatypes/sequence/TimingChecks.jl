const PULSEQ_TIME_TOL = 1e-12
const PULSEQ_DIVISION_TOL = 1e-9
const PULSEQ_ADC_DWELL_TOL = 1e-10

function _sequence_timing_from_def(def)
    missing = [key for key in PULSEQ_RASTER_DEFINITION_KEYS if !haskey(def, key)]
    isempty(missing) || error("Sequence has no Pulseq raster definitions. Use `Sequence()` or `Sequence(sys)`.")
    return _sequence_timing_from_sys(_sequence_scanner_from_def(def))
end

_sequence_timing_from_sys(sys::Scanner) = (;
    BlockDurationRaster=sys.DUR_Δt,
    GradientRasterTime=sys.GR_Δt,
    RadiofrequencyRasterTime=sys.RF_Δt,
    AdcRasterTime=sys.ADC_Δt,
    RfRingdownTime=sys.RF_ring_down_T,
    RfDeadTime=sys.RF_dead_time_T,
    AdcDeadTime=sys.ADC_dead_time_T,
)

_timing_value(timing, key, default) = hasproperty(timing, key) ? getproperty(timing, key) : default

function _check_block_fit(event_end, block_duration, block_id, label)
    event_end <= block_duration + PULSEQ_TIME_TOL && return nothing
    error("Block $block_id $label exceeds the block duration.")
end

function _check_raster_multiple(t, raster, block_id, label; tol=PULSEQ_DIVISION_TOL)
    raster > 0 || error("Raster time must be positive.")
    abs(t / raster - round(t / raster)) <= tol && return nothing
    error("Block $block_id $label ($(t) s) is not aligned to raster $(raster) s.")
end

function _check_sample_timing(start, step::Number, n, raster, block_id, label)
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

function _pulseq_uniform_interval(gr::UniformlySampledGrad)
    n = length(gr.A)
    return n > 1 ? gr.T / (n - 1) : nothing
end
function _pulseq_uniform_interval(gr::TimeShapedGrad)
    n = length(gr.A)
    lenT = length(gr.T)
    (lenT == n - 1 || lenT == n) || throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $lenT."))
    n > 1 || return nothing
    interval = gr.T[1]
    for i in 2:(n - 1)
        isapprox(gr.T[i], interval; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    end
    return interval
end

_pulseq_compact_interval(::TrapezoidalGrad, raster) = nothing
function _pulseq_compact_interval(gr::Union{UniformlySampledGrad,TimeShapedGrad}, raster)
    interval = _pulseq_uniform_interval(gr)
    isnothing(interval) && return nothing
    isapprox(gr.rise, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    isapprox(gr.fall, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    isapprox(interval, raster; rtol=0, atol=PULSEQ_TIME_TOL) && return raster
    isodd(length(gr.A)) && isapprox(interval, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) && return raster / 2
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
    compact_interval = _pulseq_compact_interval(gr, raster)
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
    compact_interval = _pulseq_compact_interval(gr, raster)
    edge_raster = isnothing(compact_interval) ? raster : raster / 2
    _check_raster_multiple(gr.delay, raster, block_id, "$name-gradient delay")
    _check_raster_multiple(gr.rise, edge_raster, block_id, "$name-gradient rise time")
    isnothing(compact_interval) && _check_sample_timing(gr.T, n, gr.rise, raster, block_id, "$name-gradient timing")
    _check_raster_multiple(gr.fall, edge_raster, block_id, "$name-gradient fall time")
    _check_block_fit(dur(gr), block_duration, block_id, "$name-gradient event")
    return nothing
end

_pulseq_compact_interval(::Union{BlockPulseRF,TimeShapedRF}, rf_raster) = nothing
function _pulseq_compact_interval(rf::UniformlySampledRF, rf_raster)
    n = length(rf.A)
    iszero(n) && return nothing
    step = n > 1 ? rf.T / (n - 1) : rf_raster
    isapprox(step, rf_raster; rtol=0, atol=PULSEQ_TIME_TOL) && return rf_raster
    isodd(n) && isapprox(step, rf_raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) && return rf_raster / 2
    return nothing
end

_pulseq_first_sample_offset(::Union{BlockPulseRF,TimeShapedRF}, rf_raster) = 0.0
function _pulseq_first_sample_offset(rf::UniformlySampledRF, rf_raster)
    return isnothing(_pulseq_compact_interval(rf, rf_raster)) ? 0.0 : rf_raster / 2
end

function _check_timing(rf::BlockPulseRF, timing, block_id, block_duration)
    raster = timing.RadiofrequencyRasterTime
    _check_raster_multiple(rf.delay, raster, block_id, "RF delay")
    _check_raster_multiple(rf.T, raster, block_id, "RF timing")
    _check_duration(rf, timing, block_id, block_duration)
    return nothing
end

function _check_timing(rf::UniformlySampledRF, timing, block_id, block_duration)
    raster = timing.RadiofrequencyRasterTime
    n = length(rf.A)
    step = n > 1 ? rf.T / (n - 1) : zero(rf.T)
    offset = _pulseq_first_sample_offset(rf, raster)
    delay = rf.delay - offset
    _check_raster_multiple(delay, raster, block_id, "RF delay")
    iszero(offset) && _check_sample_timing(0.0, step, n, raster, block_id, "RF timing")
    _check_duration(rf, timing, block_id, block_duration)
    return nothing
end

function _check_timing(rf::TimeShapedRF, timing, block_id, block_duration)
    raster = timing.RadiofrequencyRasterTime
    _check_raster_multiple(rf.delay, raster, block_id, "RF delay")
    _check_sample_timing(rf.T, length(rf.A), 0.0, raster, block_id, "RF timing")
    _check_duration(rf, timing, block_id, block_duration)
    return nothing
end

_pulseq_delay(rf, rf_raster) = rf.delay - _pulseq_first_sample_offset(rf, rf_raster)

# Compact Pulseq RF has an implicit half-raster first-sample offset.
_pulseq_duration(rf, rf_raster) = dur(rf) + _pulseq_first_sample_offset(rf, rf_raster)

function _check_duration(rf, timing, block_id, block_duration)
    rf_dead_time = _timing_value(timing, :RfDeadTime, 0.0)
    rf_ringdown_time = _timing_value(timing, :RfRingdownTime, 0.0)
    rf_end = _pulseq_duration(rf, timing.RadiofrequencyRasterTime)
    _check_block_fit(rf_end, block_duration, block_id, "RF event")
    if rf_dead_time > 0
        rf_delay = _pulseq_delay(rf, timing.RadiofrequencyRasterTime)
        rf_delay + PULSEQ_TIME_TOL >= rf_dead_time || error("Block $block_id RF delay ($(rf_delay) s) is smaller than RF dead time $(rf_dead_time) s.")
    end
    if rf_ringdown_time > 0
        _check_block_fit(rf_end + rf_ringdown_time, block_duration, block_id, "RF ring-down time")
    end
    return nothing
end

function _check_dwell(dwell, raster, block_id)
    raster > 0 || error("Raster time must be positive.")
    dwell + PULSEQ_TIME_TOL >= raster || error("Block $block_id ADC dwell ($(dwell) s) is shorter than raster $(raster) s.")
    rounded = round(dwell / raster) * raster
    abs(rounded - dwell) <= PULSEQ_ADC_DWELL_TOL && return nothing
    error("Block $block_id ADC dwell ($(dwell) s) is not aligned to raster $(raster) s.")
end

function _check_timing(adc::ADC, timing, block_id, block_duration)
    dwell = _pulseq_dwell(adc)
    _check_dwell(dwell, timing.AdcRasterTime, block_id)
    adc.delay + PULSEQ_TIME_TOL >= dwell / 2 || error("Block $block_id ADC delay ($(adc.delay) s) is smaller than dwell/2 ($(dwell / 2) s).")
    pulseq_delay = adc.delay - dwell / 2
    _check_raster_multiple(pulseq_delay, timing.RadiofrequencyRasterTime, block_id, "ADC delay")
    adc_duration = _pulseq_duration(adc)
    _check_block_fit(adc_duration, block_duration, block_id, "ADC event")
    adc_dead_time = _timing_value(timing, :AdcDeadTime, 0.0)
    if adc_dead_time > 0
        pulseq_delay + PULSEQ_TIME_TOL >= adc_dead_time || error("Block $block_id ADC delay ($(pulseq_delay) s) is smaller than ADC dead time $(adc_dead_time) s.")
        _check_block_fit(adc_duration + adc_dead_time, block_duration, block_id, "post-ADC dead time")
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
    return check_timing(seq, _sequence_timing_from_def(seq.DEF))
end

function check_timing(seq::Sequence, sys::Scanner)
    return check_timing(seq, _sequence_timing_from_sys(sys))
end

function check_timing(seq::Sequence, raster::NamedTuple)
    for i in eachindex(seq.DUR)
        _check_raster_multiple(seq.DUR[i], raster.BlockDurationRaster, i, "duration")
        rf = seq.RF[1, i]
        if is_RF_on(rf)
            _check_timing(rf, raster, i, seq.DUR[i])
        end
        axis_names = ("x", "y", "z")
        for axis in axes(seq.GR, 1)
            gr = seq.GR[axis, i]
            is_GR_on(gr) || continue
            name = axis <= length(axis_names) ? axis_names[axis] : string(axis)
            _check_timing(gr, raster, i, name, seq.DUR[i])
        end
        adc = seq.ADC[i]
        is_ADC_on(adc) && _check_timing(adc, raster, i, seq.DUR[i])
        _check_timing(seq.EXT[i], i, seq.DUR[i])
    end
    return nothing
end

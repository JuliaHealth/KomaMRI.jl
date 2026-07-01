"""
    delay(event)
    delay(event, sys)

Return the Koma event delay, or the Pulseq-facing delay when a `Scanner` is
provided. The sys-aware method accounts for RF/ADC sample-edge conventions.
"""
delay(gr::Grad) = gr.delay
delay(rf::RF) = rf.delay
delay(adc::ADC) = adc.delay
delay(::Extension) = 0.0
delay(t::Trigger) = t.delay

delay(gr::Grad, ::Scanner) = delay(gr)
delay(rf::RF, ::Scanner) = delay(rf)
function delay(rf::UniformlySampledRF, sys::Scanner)
    is_RF_on(rf) || return delay(rf)
    interval = compact_sample_interval(rf, sys.RF_Δt)
    sample_offset = interval == sys.RF_Δt ? sys.RF_Δt / 2 : dwell(rf) / 2
    delay(rf) + PULSEQ_TIME_TOL < sample_offset && return delay(rf)
    pulseq_delay = delay(rf) - sample_offset
    return interval == sys.RF_Δt ? pulseq_delay : floor_to_raster(pulseq_delay, sys.RF_Δt)
end
function delay(rf::TimeShapedRF, sys::Scanner)
    is_RF_on(rf) || return delay(rf)
    raster_delay = floor_to_raster(delay(rf), sys.RF_Δt)
    offset = delay(rf) - raster_delay
    offset > PULSEQ_TIME_TOL && return raster_delay
    sample_offset = first(rf.T) / 2
    delay(rf) + PULSEQ_TIME_TOL < sample_offset && return delay(rf)
    return floor_to_raster(delay(rf) - sample_offset, sys.RF_Δt)
end
delay(adc::ADC, sys::Scanner) = adc.delay - dwell(adc, sys) / 2
delay(ext::Extension, ::Scanner) = delay(ext)

rf_center(rf::RF, ::Scanner) = rf_center(rf)
function rf_center(rf::UniformlySampledRF, sys::Scanner)
    is_RF_on(rf) || return rf_center(rf)
    return rf_center(rf) + delay(rf) - delay(rf, sys)
end
function rf_center(rf::TimeShapedRF, sys::Scanner)
    is_RF_on(rf) || return rf_center(rf)
    return rf_center(rf) + delay(rf) - delay(rf, sys)
end

"""
    dwell(event)
    dwell(event, sys)

Return the event sample spacing. The sys-aware method is provided for API
symmetry with `dur(event, sys)` and `delay(event, sys)`.
"""
dwell(gr::Union{UniformlySampledGrad,TimeShapedGrad}, ::Scanner) = dwell(gr)
dwell(rf::Union{UniformlySampledRF,TimeShapedRF}, ::Scanner) = dwell(rf)
dwell(adc::ADC, ::Scanner) = dwell(adc)

"""
    dur(event)
    dur(event, sys)
    dur(seq, sys)

Return the Koma event duration, or the Pulseq-facing duration when a `Scanner`
is provided. RF sys-aware duration includes sample-edge timing and ring-down;
ADC sys-aware duration includes dwell-edge timing and ADC dead time.
"""
dur(gr::Grad, ::Scanner) = dur(gr)
dur(rf::RF, sys::Scanner) =
    is_RF_on(rf) ? delay(rf, sys) + sum(rf.T) + sys.RF_ring_down_time : dur(rf)
function dur(rf::UniformlySampledRF, sys::Scanner)
    is_RF_on(rf) || return dur(rf)
    shape_duration = ceil_to_raster(dur(rf) + dwell(rf) / 2 - delay(rf, sys), sys.RF_Δt)
    return delay(rf, sys) + shape_duration + sys.RF_ring_down_time
end
function dur(rf::TimeShapedRF, sys::Scanner)
    is_RF_on(rf) || return dur(rf)
    shape_duration = ceil_to_raster(dur(rf) - delay(rf, sys), sys.RF_Δt)
    return delay(rf, sys) + shape_duration + sys.RF_ring_down_time
end
dur(adc::ADC, sys::Scanner) = is_ADC_on(adc) ?
    delay(adc, sys) + adc.N * dwell(adc, sys) + sys.ADC_dead_time :
    dur(adc)
dur(ext::Extension, ::Scanner) = dur(ext)

function dur(seq::Sequence, sys::Scanner)
    duration = 0.0
    for block in eachindex(seq.DUR)
        block_duration = seq.DUR[block]
        for gr in view(seq.GR, :, block)
            block_duration = max(block_duration, dur(gr, sys))
        end
        for rf in view(seq.RF, :, block)
            block_duration = max(block_duration, dur(rf, sys))
        end
        block_duration = max(block_duration, dur(seq.ADC[block], sys))
        for ext in seq.EXT[block]
            block_duration = max(block_duration, dur(ext, sys))
        end
        duration += block_duration
    end
    return duration
end

compact_sample_interval(::TrapezoidalGrad, raster) = nothing
function compact_sample_interval(gr::UniformlySampledGrad, raster)
    interval = length(gr.A) > 1 ? dwell(gr) : nothing
    return compact_sample_interval(gr, interval, raster)
end
function compact_sample_interval(gr::TimeShapedGrad, raster)
    n = length(gr.A)
    intervals = dwell(gr)
    lenT = length(intervals)
    (lenT == n - 1 || lenT == n) || throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $lenT."))
    n > 1 || return nothing
    interval = intervals[1]
    for i in 2:(n - 1)
        isapprox(intervals[i], interval; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    end
    return compact_sample_interval(gr, interval, raster)
end

function compact_sample_interval(gr::Grad, interval, raster)
    isnothing(interval) && return nothing
    isapprox(gr.rise, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    isapprox(gr.fall, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    isapprox(interval, raster; rtol=0, atol=PULSEQ_TIME_TOL) && return raster
    isodd(length(gr.A)) && isapprox(interval, raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) && return raster / 2
    return nothing
end

compact_sample_interval(::Union{BlockPulseRF,TimeShapedRF}, rf_raster) = nothing
function compact_sample_interval(rf::UniformlySampledRF, rf_raster)
    n = length(rf.A)
    iszero(n) && return nothing
    step = n > 1 ? dwell(rf) : rf_raster
    isapprox(step, rf_raster; rtol=0, atol=PULSEQ_TIME_TOL) && return rf_raster
    isodd(n) && isapprox(step, rf_raster / 2; rtol=0, atol=PULSEQ_TIME_TOL) && return rf_raster / 2
    return nothing
end

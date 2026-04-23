"""
    A = ampls(g::Grad)
    A = ampls(r::RF)
    A = ampls(d::ADC)

Get amplitude samples of MRI sequence event.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC struct

# Returns
- `A`: (`::Vector{Number}`) vector with amplitude samples
"""
ampls(gr::TrapezoidalGrad) = is_on(gr) ? [gr.first; gr.A; gr.A; gr.last] : typeof(gr.A)[]
ampls(gr::Union{UniformlySampledGrad,TimeShapedGrad}) =
    is_on(gr) ? [gr.first; gr.A; gr.last] : eltype(gr.A)[]

function ampls(rf::RF; freq_in_phase=false)
    waveform = cis(rf.ϕ) .* rf.A
    if !is_on(rf)
        return similar(_shape_samples(waveform), 0)
    end
    A = _shape_samples(waveform)
    if freq_in_phase
        t0 = times(rf, :Δf)
        t  = times(rf)
        Interpolations.deduplicate_knots!(t0; move_knots=true)
        Interpolations.deduplicate_knots!(t; move_knots=true)
        f = linear_interpolation(t0, freqs(rf)).(t)
        phase_cycles = [0; cumtrapz(diff(t), f)]
        center_time = rf.delay + get_RF_center(rf)
        center_phase_cycles = linear_interpolation(t, phase_cycles, extrapolation_bc=Interpolations.Flat())(center_time)
        A = A .* cis.(2π .* (phase_cycles .- center_phase_cycles))
    end
	return A
end
"""
    f = freqs(r::RF)

Get frequency samples of MRI sequence event.

# Arguments
- `rf`: (`::RF`) RF struct

# Returns
- `f`: (`::Vector{Number}`) vector with frequency samples
"""
function freqs(rf::RF)
    if !is_on(rf)
        return similar(_shape_samples(rf.Δf), 0)
    end
    return _shape_samples(rf.Δf)
end
function ampls(adc::ADC)
    if !is_on(adc)
        return Bool[]
    end
    return ones(Bool, adc.N)
end


"""
    t = times(gr::Grad)
    t = times(rf::RF)
    t = times(adc::ADC)

Get time samples of MRI sequence event.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC struct

# Returns
- `t`: (`::Vector{Number}`) vector with time samples
"""

# Keep the closing knot distinct from the block boundary without changing the event shape.
function _separate_closing_knot!(t)
    length(t) >= 2 && (t[end - 1] -= MIN_RISE_TIME)
    return t
end

_gradient_times(gr::TrapezoidalGrad) = cumsum([gr.delay; gr.rise; gr.T; gr.fall])

times(gr::TrapezoidalGrad) =
    is_on(gr) ? _separate_closing_knot!(_gradient_times(gr)) : typeof(gr.delay)[]

function times(gr::UniformlySampledGrad)
    is_on(gr) || return typeof(gr.delay)[]
    return _separate_closing_knot!(_gradient_times(gr))
end

function _gradient_times(gr::UniformlySampledGrad)
    n_intervals = length(gr.A) - 1
    flat_times = n_intervals > 0 ? fill(gr.T / n_intervals, n_intervals) : typeof(gr.delay)[]
    return cumsum([gr.delay; gr.rise; flat_times; gr.fall])
end

_gradient_times(gr::TimeShapedGrad) = cumsum([gr.delay; gr.rise; gr.T; gr.fall])

times(gr::TimeShapedGrad) =
    is_on(gr) ? _separate_closing_knot!(_gradient_times(gr)) : similar(gr.T, 0)

function times(rf::RF, key::Symbol)
    if !is_on(rf)
        return typeof(rf.delay)[]
    end
    ts0 = key === :A ? _shape_times(rf.A, rf.T) :
        key === :Δf ? _shape_times(rf.Δf, rf.T) :
        throw(ArgumentError("Unsupported RF key $key"))
    ts = rf.delay .+ ts0
    t = [rf.delay; ts; ts[end]]
    return _separate_closing_knot!(t)
end
times(rf::RF) = times(rf, :A)
function times(adc::ADC)
    if !is_on(adc)
        return range(0.0,0.0,0) #empty
    end
	if adc.N != 1
        t = range(0.0, adc.T; length=adc.N) .+ adc.delay
    else
        t = range(adc.T/2, adc.T/2, 1) .+ adc.delay
    end
    return t
end

is_on(x::Number) = !iszero(x)
is_on(x::AbstractArray{<:Number}) = any(!iszero, x)
is_on(x::Grad) = is_on(x.A)
is_on(x::RF)   = is_on(x.A)
is_on(x::ADC)  = x.N > 0

is_GR_on(x::Grad) = is_on(x)
is_RF_on(x::RF) = is_on(x)
is_ADC_on(x::ADC) = is_on(x)

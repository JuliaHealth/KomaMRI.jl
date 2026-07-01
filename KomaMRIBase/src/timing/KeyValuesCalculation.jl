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

function area(gr::TrapezoidalGrad)
    is_on(gr) || return 0.0
    return gr.A * gr.T +
        (gr.first + gr.A) * gr.rise / 2 +
        (gr.A + gr.last) * gr.fall / 2
end

function area(gr::Union{UniformlySampledGrad,TimeShapedGrad})
    A = ampls(gr)
    isempty(A) && return 0.0
    return trapz(diff(times(gr; separate_closing_knot=false)), A)
end

function area(rf::RF)
    A = ampls(rf)
    isempty(A) && return zero(eltype(A))
    return trapz(diff(times(rf; separate_closing_knot=false)), A)
end

"""
    Δt = dwell(event)

Return the event sample spacing in seconds. Time-shaped gradient and RF events
return their stored interval vector.
"""
dwell(gr::UniformlySampledGrad) = length(gr.A) <= 1 ? sum(gr.T) : gr.T / (length(gr.A) - 1)
dwell(gr::TimeShapedGrad) = gr.T
dwell(rf::UniformlySampledRF) = length(rf.A) <= 1 ? sum(rf.T) : rf.T / (length(rf.A) - 1)
dwell(rf::TimeShapedRF) = rf.T
dwell(adc::ADC) = adc.N <= 1 ? adc.T : adc.T / (adc.N - 1)

function _with_frequency_phase(rf, A, freq_in_phase)
    if !is_on(rf)
        return similar(A, 0)
    end
    if freq_in_phase
        t0 = freq_times(rf)
        t  = times(rf)
        Interpolations.deduplicate_knots!(t0; move_knots=true)
        Interpolations.deduplicate_knots!(t; move_knots=true)
        f = linear_interpolation(t0, freqs(rf)).(t)
        phase_cycles = [0; cumtrapz(diff(t), f)]
        center_time = rf.delay + rf_center(rf)
        center_phase_cycles = linear_interpolation(t, phase_cycles, extrapolation_bc=Interpolations.Flat())(center_time)
        A = A .* cis.(2π .* (phase_cycles .- center_phase_cycles))
    end
	return A
end

function ampls(rf::BlockPulseRF; freq_in_phase=false)
    waveform = cis(rf.ϕ) * rf.A
    A = [zero(waveform), waveform, waveform, zero(waveform)]
    return _with_frequency_phase(rf, A, freq_in_phase)
end
function ampls(rf::Union{UniformlySampledRF,TimeShapedRF}; freq_in_phase=false)
    waveform = cis(rf.ϕ) .* rf.A
    A = [zero(eltype(waveform)); waveform; zero(eltype(waveform))]
    return _with_frequency_phase(rf, A, freq_in_phase)
end
"""
    f = freqs(r::RF)

Get frequency samples of MRI sequence event.

# Arguments
- `rf`: (`::RF`) RF struct

# Returns
- `f`: (`::Vector{Number}`) vector with frequency samples
"""
function freqs(rf::FrequencyModulatedRF)
    f = [zero(eltype(rf.Δf)); rf.Δf; zero(eltype(rf.Δf))]
    if !is_on(rf)
        return similar(f, 0)
    end
    return f
end
function freqs(rf::RF)
    f = [zero(rf.Δf), rf.Δf, rf.Δf, zero(rf.Δf)]
    if !is_on(rf)
        return similar(f, 0)
    end
    return f
end

"""
    ψ = rf_frame_phase(rf::RF)

RF rotating-frame phase samples from frequency modulation, with `ψ` zeroed at the RF center.
"""
function rf_frame_phase(rf::RF)
    f = freqs(rf)
    isempty(f) && return Float64[]
    t = freq_times(rf)
    Interpolations.deduplicate_knots!(t; move_knots=true)
    center_time = rf.delay + rf_center(rf)
    t_aug = sort!([t; center_time])
    f_aug = linear_interpolation(t, f, extrapolation_bc=Interpolations.Flat()).(t_aug)
    ψ_aug = -2π .* [0.0; cumtrapz(diff(t_aug), f_aug)]
    ψ = ψ_aug[searchsortedfirst.(Ref(t_aug), t)]
    ψ_center = ψ_aug[searchsortedfirst(t_aug, center_time)]
    ψ .-= ψ_center
    ψ[1] = zero(eltype(ψ))
    ψ[end] = zero(eltype(ψ))
    return ψ
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

function _reseparate_closing_knot!(t)
    length(t) >= 2 || return t
    t[end - 1] = min(t[end - 1], t[end] - max(MIN_RISE_TIME, eps(t[end])))
    return t
end

function times(gr::TrapezoidalGrad; separate_closing_knot=true)
    is_on(gr) || return typeof(gr.delay)[]
    t = cumsum([gr.delay; gr.rise; gr.T; gr.fall])
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function times(gr::UniformlySampledGrad; separate_closing_knot=true)
    is_on(gr) || return typeof(gr.delay)[]
    n_intervals = length(gr.A) - 1
    flat_times = n_intervals > 0 ? fill(gr.T / n_intervals, n_intervals) : typeof(gr.delay)[]
    t = cumsum([gr.delay; gr.rise; flat_times; gr.fall])
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function times(gr::TimeShapedGrad; separate_closing_knot=true)
    is_on(gr) || return similar(gr.T, 0)
    t = cumsum([gr.delay; gr.rise; gr.T; gr.fall])
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function times(rf::BlockPulseRF; separate_closing_knot=true)
    is_on(rf) || return typeof(rf.delay)[]
    rf_duration = sum(rf.T)
    t = rf.delay .+ [zero(rf_duration), zero(rf_duration), rf_duration, rf_duration]
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function times(rf::UniformlySampledRF; separate_closing_knot=true)
    is_on(rf) || return typeof(rf.delay)[]
    ts = rf.delay .+ range(zero(rf.T), rf.T; length=length(rf.A))
    t = [rf.delay; ts; ts[end]]
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function times(rf::TimeShapedRF; separate_closing_knot=true)
    n = length(rf.A)
    n == 0 && return typeof(rf.delay)[]
    ts0 = length(rf.T) == n - 1 ? cumsum([zero(eltype(rf.T)); rf.T]) :
        length(rf.T) == n ? cumsum([zero(eltype(rf.T)); rf.T[1:(end - 1)]]) :
        throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $(length(rf.T))."))
    ts = rf.delay .+ ts0
    t = [rf.delay; ts; ts[end]]
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

"""
    t = freq_times(rf::RF)

Get frequency-modulation sample times of an RF event.
"""
function freq_times(rf::RF; separate_closing_knot=true)
    is_on(rf) || return typeof(rf.delay)[]
    rf_duration = sum(rf.T)
    t = rf.delay .+ [zero(rf_duration), zero(rf_duration), rf_duration, rf_duration]
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function freq_times(rf::UniformlySampledFrequencyModulatedRF; separate_closing_knot=true)
    is_on(rf) || return typeof(rf.delay)[]
    ts = rf.delay .+ range(zero(rf.T), rf.T; length=length(rf.Δf))
    t = [rf.delay; ts; ts[end]]
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function freq_times(rf::TimeShapedFrequencyModulatedRF; separate_closing_knot=true)
    n = length(rf.Δf)
    n == 0 && return typeof(rf.delay)[]
    ts0 = length(rf.T) == n - 1 ? cumsum([zero(eltype(rf.T)); rf.T]) :
        length(rf.T) == n ? cumsum([zero(eltype(rf.T)); rf.T[1:(end - 1)]]) :
        throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $(length(rf.T))."))
    ts = rf.delay .+ ts0
    t = [rf.delay; ts; ts[end]]
    return separate_closing_knot ? _separate_closing_knot!(t) : t
end

function times(adc::ADC)
    if !is_on(adc)
        return range(0.0,0.0,0) #empty
    end
	if adc.N != 1
        t = range(0.0, adc.T; length=adc.N) .+ adc.delay
    else
        t = range(0.0, 0.0, 1) .+ adc.delay
    end
    return t
end

is_on(x::Number) = !iszero(x)
is_on(x::AbstractArray{<:Number}) = any(!iszero, x)
is_on(x::Grad) = is_on(x.A) || is_on(x.first) || is_on(x.last)
is_on(x::RF)   = is_on(x.A)
is_on(x::ADC)  = x.N > 0

is_GR_on(x::Grad) = is_on(x)
is_RF_on(x::RF) = is_on(x)
is_ADC_on(x::ADC) = is_on(x)

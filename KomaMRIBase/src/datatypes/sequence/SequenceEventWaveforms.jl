# ===========================================================================
# 1. Sequence Event Waveforms
# ===========================================================================
#
# MRI sequence events (RF, gradients, ADC) are converted here to local event
# waveforms: time samples plus RF amplitude, gradient amplitude, frequency, or
# acquisition gates.

# -- 1.1. Event activity predicates -----------------------------------------
is_on(x::Number) = !iszero(x)
is_on(x::AbstractArray{<:Number}) = any(!iszero, x)
is_on(x::Grad) = is_on(x.A)
is_on(x::RF)   = is_on(x.A)
is_on(x::ADC)  = x.N > 0

is_GR_on(x::Grad) = is_on(x)
is_RF_on(x::RF) = is_on(x)
is_ADC_on(x::ADC) = is_on(x)

# -- 1.2. Event amplitudes and RF frequency offsets --------------------------
"""
    A = ampls(g::Grad)
    A = ampls(r::RF)
    A = ampls(d::ADC)

Get amplitude samples of MRI sequence event.
"""
ampls(gr::TrapezoidalGrad) = is_on(gr) ? [gr.first; gr.A; gr.A; gr.last] : typeof(gr.A)[]
function ampls(gr::Union{UniformlySampledGrad,TimeShapedGrad})
    is_on(gr) || return eltype(gr.A)[]
    A = length(gr.A) == 1 ? [only(gr.A), only(gr.A)] : gr.A
    return [gr.first; A; gr.last]
end

rf_center_time(rf::RF) = rf.delay + rf_center(rf)

function ampls(rf::BlockPulseRF; freq_in_phase=false)
    A0 = cis(rf.ϕ) * rf.A
    is_on(rf) || return typeof(A0)[]
    A = [zero(A0), A0, A0, zero(A0)]
    if freq_in_phase
        t  = times(rf)
        Δf = (t=times(rf, :Δf)[2:(end - 1)], A=freqs(rf)[2:(end - 1)])
        A .*= cis.(rf_phase_at_times(Δf, t, rf_center_time(rf)))
    end
    return A
end

function ampls(rf::RF; freq_in_phase=false)
    A = collect(cis(rf.ϕ) .* rf.A)
    is_on(rf) || return similar(A, 0)
    length(A) == 1 && (A = [only(A), only(A)])
    A = [zero(eltype(A)); A; zero(eltype(A))]
    if freq_in_phase
        t  = times(rf)
        Δf = (t=times(rf, :Δf)[2:(end - 1)], A=freqs(rf)[2:(end - 1)])
        A .*= cis.(rf_phase_at_times(Δf, t, rf_center_time(rf)))
    end
    return A
end

"""
    f = freqs(r::RF)

Get frequency samples of MRI RF event.
"""
function freqs(rf::RF)
    is_on(rf) || return typeof(rf.Δf)[]
    Δf = fill(rf.Δf, length(freq_times(rf)))
    Δf[begin] = zero(rf.Δf)
    Δf[end] = zero(rf.Δf)
    return Δf
end

function freqs(rf::FrequencyModulatedRF)
    is_on(rf) || return eltype(rf.Δf)[]
    Δf = collect(rf.Δf)
    length(Δf) == 1 && (Δf = [only(Δf), only(Δf)])
    return [zero(eltype(Δf)); Δf; zero(eltype(Δf))]
end

function ampls(adc::ADC)
    if !is_on(adc)
        return Bool[]
    end
    return ones(Bool, adc.N)
end

# -- 1.3. Local event timing -------------------------------------------------
"""
    t = times(gr::Grad)
    t = times(rf::RF)
    t = times(adc::ADC)

Get time samples of MRI sequence event.
"""
times(gr::TrapezoidalGrad) =
    is_on(gr) ? cumsum([gr.delay; gr.rise; gr.T; gr.fall]) : typeof(gr.delay)[]

function times(gr::UniformlySampledGrad)
    is_on(gr) || return typeof(gr.delay)[]
    n_intervals = length(gr.A) - 1
    flat_times = n_intervals > 0 ? fill(gr.T / n_intervals, n_intervals) : [gr.T]
    return cumsum([gr.delay; gr.rise; flat_times; gr.fall])
end

times(gr::TimeShapedGrad) =
    is_on(gr) ? cumsum([gr.delay; gr.rise; gr.T; gr.fall]) : similar(gr.T, 0)

function times(rf::BlockPulseRF)
    is_on(rf) || return typeof(rf.delay)[]
    return [rf.delay, rf.delay, dur(rf), dur(rf)]
end

function times(rf::UniformlySampledRF)
    is_on(rf) || return typeof(rf.delay)[]
    t = length(rf.A) == 1 ? [rf.delay, dur(rf)] : collect(range(rf.delay, dur(rf); length=length(rf.A)))
    return [first(t); t; last(t)]
end

function times(rf::TimeShapedRF)
    is_on(rf) || return typeof(rf.delay)[]
    t = rf.delay .+ cumsum(rf.T)
    t = length(rf.T) == length(rf.A) - 1 ? [rf.delay; t] : t
    return [first(t); t; last(t)]
end

"""
    t = freq_times(rf::RF)

Get RF frequency modulation time samples.
"""
function freq_times(rf::RF)
    is_on(rf) || return typeof(rf.delay)[]
    return [rf.delay, rf.delay, dur(rf), dur(rf)]
end

function freq_times(rf::UniformlySampledFrequencyModulatedRF)
    is_on(rf) || return typeof(rf.delay)[]
    t = length(rf.Δf) == 1 ? [rf.delay, dur(rf)] : collect(range(rf.delay, dur(rf); length=length(rf.Δf)))
    return [first(t); t; last(t)]
end

function freq_times(rf::TimeShapedFrequencyModulatedRF)
    is_on(rf) || return typeof(rf.delay)[]
    t = rf.delay .+ cumsum(rf.T)
    t = length(rf.T) == length(rf.Δf) - 1 ? [rf.delay; t] : t
    return [first(t); t; last(t)]
end

function times(rf::RF, key::Symbol)
    key === :A && return times(rf)
    key === :Δf && return freq_times(rf)
    throw(ArgumentError("Unsupported RF key $key"))
end
times(rf::RF, key) = throw(ArgumentError("Unsupported RF key $key"))

function times(adc::ADC)
    if !is_on(adc)
        return Float64[]
    end
    adc.N == 1 && return [adc.delay + adc.T / 2]
    return collect(range(0.0, adc.T; length=adc.N) .+ adc.delay)
end

# -- 1.4. RF rotating-frame phase -------------------------------------------
function rf_phase_at_times(Δf, t, center_time)
    ψ = zeros(length(t))
    (isempty(Δf.t) || isempty(t) || isnothing(center_time)) && return ψ
    lo = searchsortedfirst(t, first(Δf.t))
    hi = searchsortedlast(t, last(Δf.t))
    lo <= hi || return ψ
    idx = lo:hi
    integral = trapezoidal_antiderivative(Δf)
    center_integral = waveform_integral_at(Δf, integral, center_time)
    for i in idx
        ψ[i] = -2π * (waveform_integral_at(Δf, integral, t[i]) - center_integral)
    end
    return ψ
end

function trapezoidal_antiderivative(samples)
    integral = zeros(length(samples.t))
    for i in 1:(length(samples.t) - 1)
        integral[i + 1] = integral[i] + (samples.t[i + 1] - samples.t[i]) * (samples.A[i] + samples.A[i + 1]) / 2
    end
    return integral
end

function waveform_integral_at(samples, integral, t)
    i = searchsortedlast(samples.t, t)
    i <= 0 && return zero(eltype(integral))
    i >= length(samples.t) && return integral[end]
    t0, t1 = samples.t[i], samples.t[i + 1]
    t1 == t0 && return integral[i]
    At = samples.A[i] + (samples.A[i + 1] - samples.A[i]) * ((t - t0) / (t1 - t0))
    return integral[i] + (t - t0) * (samples.A[i] + At) / 2
end

# -- 1.5. Event waveform sample pairs ---------------------------------------
event_samples(gr::Grad) = (t=times(gr), A=ampls(gr))
event_samples(rf::RF; freq_in_phase=false) = (t=times(rf), A=ampls(rf; freq_in_phase))
event_samples(rf::RF, ::Val{:A}; freq_in_phase=false) = event_samples(rf; freq_in_phase)
event_samples(rf::RF, ::Val{:Δf}) = (t=times(rf, :Δf), A=freqs(rf))

function event_samples(rf::RF, ::Val{:ψ})
    Δf = event_samples(rf, Val(:Δf))
    ψ = isempty(Δf.t) ? similar(Δf.A, 0) : rf_phase_at_times(Δf, Δf.t, rf_center_time(rf))
    return (t=Δf.t, A=ψ)
end

event_samples(rf::RF, key::Symbol; kwargs...) = event_samples(rf, Val(key); kwargs...)
event_samples(rf::RF, key; kwargs...) = throw(ArgumentError("Unsupported RF key $key"))
event_samples(adc::ADC) = (t=times(adc), A=ampls(adc))

"""
    adc = ADC(N, T)
    adc = ADC(N, T, delay)
    adc = ADC(N, T, delay, Δf, ϕ)

The ADC struct represents the Analog to Digital Converter (ADC) of a sequence event.

# Arguments
- `N`: (`::Int64`) number of acquired samples
- `T`: (`::Float64`, [`s`]) duration to acquire the samples
- `delay`: (`::Float64`, [`s`]) delay time to start the acquisition
- `Δf`: (`::Float64`, [`Hz`]) delta frequency. It is meant to compensate RF pulse phases
- `ϕ`: (`::Float64`, `[rad]`) phase. It is meant to compensate RF pulse phases

# Returns
- `adc`: (`::ADC`) ADC struct

# Examples
```julia-repl
julia> adc = ADC(16, 1, 0.1)

julia> seq = Sequence(); seq += adc; plot_seq(seq)
```
"""
mutable struct ADC
    N::Int64
    T::Float64
    delay::Float64
    Δf::Float64
    ϕ::Float64
    function ADC(N, T, delay, Δf, ϕ)
        T < 0 || delay < 0 ? error("ADC timings must be positive.") : new(N, T, delay, Δf, ϕ)
    end
    function ADC(N, T, delay)
		T < 0 || delay < 0 ? error("ADC timings must be positive.") : new(N, T, delay, 0, 0)
    end
    function ADC(N, T)
		T < 0 ? error("ADC timings must be positive.") : new(N, T, 0, 0, 0)
    end
end

# ADC comparison
Base.isapprox(adc1::ADC, adc2::ADC) = begin
    return all(length(getfield(adc1, k)) ≈ length(getfield(adc2, k)) for k ∈ fieldnames(ADC))
        all(getfield(adc1, k) ≈ getfield(adc2, k) for k ∈ fieldnames(ADC))
end

"""
    y = getproperty(x::Vector{ADC}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the ADC vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{ADC}`) vector of ADC structs
- `f`: (`::Symbol`, opts: [`:N`, `:T`, `:delay`, `:Δf`, `:ϕ`, `:dur`]) input symbol that
    represents a property of the ADC structs

# Returns
- `y`: (`::Vector{Any}`) vector with the property defined by the `f` for all elements of
    the ADC vector `x`
"""
getproperty(x::Vector{ADC}, f::Symbol) = begin
    if f == :dur
		dur.(x)
    elseif f in fieldnames(ADC)
        getfield.(x, f)
    else
        getfield(x, f)
    end
end

"""
    times = get_adc_sampling_times(seq)

Returns an array of times when the samples of the sequence `seq` are acquired.

# Arguments
- `seq`: (`::Sequence`) sequence struct

# Returns
- `times`: (`::Vector{Float64}`, `[s]`) time array when samples are acquired
"""
function get_adc_sampling_times(seq)
    T0 = get_block_start_times(seq)
    t = zeros(Float64, sum(seq.ADC.N))
    idx = 1
    for i = 1:length(seq)
        if is_ADC_on(seq[i])
            N = seq.ADC[i].N
            t[idx:idx+N-1] .= times(seq.ADC[i]) .+ T0[i]
            idx += N
        end
    end
    return t
end

"""
    comp = get_adc_phase_compensation(seq)

Returns an array of phase compensation factors, ``\\exp(-\\mathrm{i}\\varphi)``, which are
used to compensate the acquired signal ``S`` by applying the operation
``S_{\\mathrm{comp}} = S \\exp(-\\mathrm{i}\\varphi)`` after the simulation. This compensation
is necessary because the signal typically exhibits a phase offset of ``\\varphi`` following
RF excitation with a phase of ``\\varphi``. Such pulses are commonly employed in sequences
involving RF spoiling.

# Arguments
- `seq`: (`::Sequence`) sequence struct

# Returns
- `comp`: (`::Vector{Complex}`, `[rad]`) array of phase compensations for every acquired sample
"""
function get_adc_phase_compensation(seq)
    phase = ComplexF32[]
    for i in 1:length(seq)
        if is_ADC_on(seq[i])
            N = seq.ADC[i].N
            ϕ = seq.ADC[i].ϕ
            aux = ones(N) .* exp(-1im * ϕ)
            append!(phase, aux)
        end
    end
    return phase
end

dur(adc::ADC) = adc.delay + adc.T

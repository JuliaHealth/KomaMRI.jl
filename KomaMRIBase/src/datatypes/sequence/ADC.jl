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
mutable struct ADC <: MRISequenceEvent
    N::Int64
    T::Float64
    delay::Float64
    Δf::Float64
    ϕ::Float64
    function ADC(N, T, delay, Δf, ϕ)
        return if T < 0 || delay < 0
            error("ADC timings must be positive.")
        else
            new(N, T, delay, Δf, ϕ)
        end
    end
    function ADC(N, T, delay)
        return if T < 0 || delay < 0
            error("ADC timings must be positive.")
        else
            new(N, T, delay, 0, 0)
        end
    end
    function ADC(N, T)
        return T < 0 ? error("ADC timings must be positive.") : new(N, T, 0, 0, 0)
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
    times = Float64[]
    for i in 1:length(seq)
        if is_ADC_on(seq[i])
            t = time(seq.ADC[i]) .+ T0[i]
            append!(times, t)
        end
    end
    return times
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
dur(adc::Vector{ADC}) = dur.(adc)

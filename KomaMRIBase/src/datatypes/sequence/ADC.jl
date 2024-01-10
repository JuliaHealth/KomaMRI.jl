"""
    adc = ADC(N, T)
    adc = ADC(N, T, delay)
    adc = ADC(N, T, delay, Δf, ϕ)

The ADC struct represents the Analog to Digital Converter (ADC) event of a sequence.

# Arguments
- `N`: (`::Integer`) number of acquired samples
- `T`: (`::Real`, [`s`]) duration to acquire the samples
- `delay`: (`::Real`, [`s`]) delay time to start the acquisition
- `Δf`: (`::Real`, [`Hz`]) delta frequency. It is meant to compensate RF pulse phases
- `ϕ`: (`::Real`, `[rad]`) phase. It is meant to compensate RF pulse phases

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

Overloads Base.getproperty() to facilitate direct access to properties of the ADC vector `x`
without the need for elementwise iteration.

# Arguments
- `x`: (`::Vector{ADC}`) vector of ADC structs
- `f`: (`::Symbol`, opts: [`:N`, `:T`, `:delay`, `:Δf`, `:ϕ`, `:dur`]) input symbol that
    represents a property of the ADC structs

# Returns
- `y`: (`::Vector{Any}`) vector with the property defined by `f` for all elements of the ADC
    vector `x`
"""
getproperty(x::Vector{ADC}, f::Symbol) = begin
    if f == :dur
		T, delay = x.T, x.delay
		ΔT = T .+ delay
		ΔT
    else
        getproperty.(x, f)
    end
end

"""
    times = get_adc_sampling_times(seq)

Returns an array of times when the samples of the sequence `seq` are acquired.

# Arguments
- `seq`: (`::Sequence`) sequence struct

# Returns
- `times`: (`::Vector{Float64}`, `[s]`) time array during sample acquisiton
"""
function get_adc_sampling_times(seq)
    T0 = get_block_start_times(seq)
    times = Float64[]
    for i = 1:length(seq)
        if is_ADC_on(seq[i])
            δ = seq.ADC[i].delay
            T = seq.ADC[i].T
            N = seq.ADC[i].N
            if N != 1
                t = range(0, T; length=N).+T0[i].+δ #range(0,T,N) works in Julia 1.7
            else
                t = [T/2].+T0[i].+δ #range(0,T,N) works in Julia 1.7
            end
            append!(times, t)
        end
    end
    return times
end

"""
    phase = get_adc_phase_compensation(seq)

Returns the array of phases for each acquired sample in the sequence `seq`. This function is
particularly useful for compensating the phase when the RF pulse also has a phase. Refer to
the end of the [`run_sim_time_iter`](@ref) function to see its usage.

# Arguments
- `seq`: (`::Sequence`) sequence struct

# Returns
- `phase`: (`::Vector{Complex{Float64}}`, `[rad]`) array of phases for each acquired sample
"""
function get_adc_phase_compensation(seq)
  phase = ComplexF32[]
  for i = 1:length(seq)
      if is_ADC_on(seq[i])
          N = seq.ADC[i].N
          ϕ = seq.ADC[i].ϕ
          aux = ones(N) .* exp(-1im*ϕ)
          append!(phase, aux)
      end
  end
  return phase
end

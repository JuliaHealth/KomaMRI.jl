"""
    adc = ADC(N, T)
    adc = ADC(N, T, delay)
    adc = ADC(N, T, delay, Δf, ϕ)

The ADC struct.

# Arguments
- `N`: (`::Int64`) number of acquired samples
- `T`: (`::Float64`, [`s`]) duration to acquire the samples
- `delay`: (`::Float64`, [`s`]) delay time to start the acquisition
- `Δf`: (`::Float64`, [`Hz`]) delta frequency. It's meant to compensate RF pulse phases.
    It is used internally by the [`read_ADC`](@ref) function
- `ϕ`: (`::Float64`, `[rad]`) phase. It's meant to compensate RF pulse phases. It is
    used internally by the [`read_ADC`](@ref) function

# Returns
- `adc`: (`::ADC`) ADC struct
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

"""
    y = getproperty(x::Vector{ADC}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the ADC vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{ADC}`) vector of ADC structs
- `f`: (`::Symbol`, opts: [`:N`, `:T`, `:delay`, `:Δf`, `:ϕ`, `:dur`]) input symbol that
    represents a property of the ACD structs

# Returns
- `y`: (`::Vector{Any}`) vector with the property defined by the `f` for all elements of
    the ADC vector `x`
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
- `times`: (`::Vector{Float64}`, `[s]`) time array when samples are acquired
"""
function get_adc_sampling_times(seq)
    T0 = cumsum([0; durs(seq)], dims=1)
    times = Float64[]
    for i = 1:length(seq)
        if is_ADC_on(seq[i])
            δ = seq.ADC[i].delay
            T = seq.ADC[i].T
            N = seq.ADC[i].N
            t = range(0, T; length=N).+T0[i].+δ #range(0,T,N) works in Julia 1.7
            append!(times, t)
        end
    end
    return times
end

"""
    phase = get_adc_phase_compensation(seq)

Returns the array of phases for every acquired sample in the sequence `seq`.

!!! note
    This function is useful to compensate the phase when the RF pulse has a phase too. Refer
    to the end of the [`run_sim_time_iter`](@ref) function to see its usage.

# Arguments
- `seq`: (`::Sequence`) sequence struct

# Returns
- `phase`: (`::Vector{Complex{Int64}}`, `[rad]`) array of phases for every acquired sample
"""
function get_adc_phase_compensation(seq)
  phase = ComplexF64[]
  for i = 1:length(seq)
      if is_ADC_on(seq[i])
          N = seq.ADC[i].N
          ϕ = seq.ADC[i].ϕ
          aux = ones(N) .* exp(1im*ϕ)
          append!(phase, aux)
      end
  end
  return phase
end

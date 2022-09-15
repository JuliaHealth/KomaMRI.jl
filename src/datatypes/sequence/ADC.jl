"""
    adc = ADC(N, T)
    adc = ADC(N, T, delay)
    adc = ADC(N, T, delay, Δf, ϕ)

The ADC struct.

# Arguments
- `N`: (`::Int64`) the number of acquired samples
- `T`: (`::Float64`, [`s`]) the duration to acquire the samples
- `delay`: (`::Float64`, [`s`]) the delay time to start the acquisition
- `Δf`: (`::Float64`, [`Hz`]) the delta frequency. It's meant to compensate RF pulse phases.
    It is used internally by the [`read_ADC`](@ref) function
- `ϕ`: (`::Float64`, `[rad]`) the phase. It's meant to compensate RF pulse phases. It is
    used internally by the [`read_ADC`](@ref) function

# Returns
- `adc`: (`::ADC`) the ADC struct

# Examples
```julia-repl
julia> d1, d2, d3 = 0.8, 0.4, 0.8;

julia> matrixGrads = [Grad(0, d1) Grad(0, d2) Grad(0, d3)];

julia> matrixRFs = [RF(0, d1) RF(0, d2) RF(0, d3)];

julia> vectorADCs = [ADC(0, d1); ADC(0, d2); ADC(9, d3)];

julia> seq = Sequence(matrixGrads, matrixRFs, vectorADCs)
Sequence[ τ = 2000.0 ms | blocks: 3 | ADC: 1 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
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

"""
    y = getproperty(x::Vector{ADC}, f::Symbol)

Overchages Base.getproperty(). It is meant to access properties of the ADC vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{ADC}`) the vector of ADC structs
- `f`: (`::Symbol`, opts: [`:N`, `:T`, `:delay`, `:Δf`, `:ϕ`, `:dur`]) the input symbol that
    represents a property of the ACD structs

# Returns
- `y`: (`::Vector{Any}`) the vector with the property defined by the `f` for all elements of
    the ADC vector `x`

``` julia-repl
julia> ADCs = [ADC(16, 8, 2); ADC(8, 4, 6); ADC(4, 2, 8)]
3-element Vector{ADC}:
 ADC(16, 8.0, 2.0, 0.0, 0.0)
 ADC(8, 4.0, 6.0, 0.0, 0.0)
 ADC(4, 2.0, 8.0, 0.0, 0.0)

julia> getproperty(ADCs, :dur)
3-element Vector{Float64}:
 10.0
 10.0
 10.0
```
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
    times = get_sample_times(seq)

Returns an array of times when the samples of the sequence `seq` are acquired.

# Arguments
- `seq`: (`::Sequence`) the sequence struct

# Returns
- `times`: (`::Vector{Float64}`, `[s]`) the time array when samples are acquired
"""
function get_sample_times(seq)
    T0 = cumsum([0; durs(seq)], dims=1)
    times = []
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
    phase = get_sample_phase_compensation(seq)

Returns the array of phases for every acquired sample in the sequence `seq`.

!!! note
    This function is useful to compensate the phase when the RF pulse has a phase too. Refer
    to the end of the [`run_sim_time_iter`](@ref) function to see its usage.

# Arguments
- `seq`: (`::Sequence`) the sequence struct

# Returns
- `phase`: (`::Vector{Complex{Int64}}`, `[rad]`) the array of phases for every acquired
    sample
"""
function get_sample_phase_compensation(seq)
  phase = []
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

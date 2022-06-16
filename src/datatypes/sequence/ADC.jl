mutable struct ADC
    N::Int64
    T::Float64
    delay::Float64
    Δf::Float64
    ϕ::Float64
    function ADC(A,T,delay,Δf,ϕ)
      T < 0 || delay < 0 ? error("ADC timings must be positive.") : new(A, T, delay, Δf, ϕ)
    end
    function ADC(A,T,delay)
		T < 0 || delay < 0 ? error("ADC timings must be positive.") : new(A, T, delay, 0, 0)
    end
    function ADC(A,T)
		T < 0 ? error("ADC timings must be positive.") : new(A, T, 0, 0, 0)
    end
end

getproperty(x::Vector{ADC}, f::Symbol) = begin
  if f == :dur
		T, delay = x.T, x.delay
		ΔT = T .+ delay
		ΔT
  else
    getproperty.(x,f)
  end
end

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
    times
end

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
  phase
end
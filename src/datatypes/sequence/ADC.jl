mutable struct ADC
    N::Int64
    T::Float64
    delay::Float64
    #freq?
    #phase?
    function ADC(A,T,delay)
		T < 0 || delay < 0 ? error("ADC timings must be positive.") : new(A, T, delay)
    end
    function ADC(A,T)
		T < 0 ? error("ADC timings must be positive.") : new(A, T, 0)
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
            t = range(T0[i].+δ, T0[i].+T.+δ, length=N)
            append!(times, t)
        end
    end
    times
end
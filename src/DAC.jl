mutable struct DAC
    N::Int64
    T::Float64
end
getproperty(x::Vector{DAC}, f::Symbol) = getproperty.(x,f)

function get_sample_times(seq)
    T = seq.DAC.T
    t0 = cumsum([0; T[1:end-1]])
    tf = t0 .+ T
    N = seq.DAC.N
    times = []
    for i = 1:length(T)
        if N[i] > 0
            t = range(t0[i], tf[i], length=N[i])
            append!(times, t)
        end
    end
    times
end
"""
    rf = RF(A, T)
    rf = RF(A, T, Δf)
    rf = RF(A, T, Δf, delay)

The RF struct represents a Radio Frequency excitation event of a sequence.

# Arguments
- `A`: (`::Complex`, `[T]`) RF complex amplitud modulation (AM), ``B_1(t) = |B_1(t)|
    e^{i\\phi(t)} = B_{1}(t) + iB_{1,y}(t) ``
- `T`: (`::Real`, [`s`]) RF duration
- `Δf`: (`::Real` or `::Vector`, [`Hz`]) RF frequency difference with respect to the Larmor
    frequency. This can be a number but also a vector to represent frequency modulated
    signals (FM).
- `delay`: (`::Real`, [`s`]) RF delay time

# Returns
- `rf`: (`::RF`) RF struct

# Examples
```julia-repl
julia> rf = RF(1, 1, 0, 0.2)

julia> seq = Sequence(); seq += rf; plot_seq(seq)
```
"""
mutable struct RF
	A
	T
	Δf
	delay::Real
	function RF(A, T, Δf, delay)
        any(T .< 0) || delay < 0 ? error("RF timings must be non-negative.") : new(A, T, Δf, delay)
    end
	function RF(A, T, Δf)
        any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, Δf, 0.)
    end
	function RF(A, T)
        any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, 0., 0.)
    end
end

# Display on the REPL
Base.show(io::IO, x::RF) = begin
	r(x) = round.(x,digits=4)
	compact = get(io, :compact, false)
	if !compact
		wave = length(x.A) == 1 ? r(x.A*1e6) : "∿"
		print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"RF($(wave) uT, $(r(sum(x.T)*1e3)) ms, $(r(x.Δf)) Hz)")
	else
		wave = length(x.A) == 1 ? "⊓" : "∿"
		print(io, (sum(abs.(x.A)) > 0 ? wave : "⇿")*"($(r((x.delay+sum(x.T))*1e3)) ms)")
	end
end

"""
    y = getproperty(x::Vector{RF}, f::Symbol)
    y = getproperty(x::Matrix{RF}, f::Symbol)

Overloads Base.getproperty(). This function is designed to directly access properties of the
RF array x without the need for elementwise iteration.

# Arguments
- `x`: (`::Vector{RF}` or `::Matrix{RF}`) vector or matrix of RF structs
- `f`: (`::Symbol`, opts: [`:A`, `:Bx`, `:By`, `:T`, `:Δf`, `:delay` and `:dur`]) input
    symbol representing a property of the vector or matrix of RF structs

# Returns
- `y`: (`::Vector{Any}` or `::Matrix{Any}`) vector containing the property defined by the
    symbol `f` for all elements of the RF vector or matrix `x`
"""
getproperty(x::Vector{RF}, f::Symbol) = getproperty.(x,f)
getproperty(x::Matrix{RF}, f::Symbol) = begin
	if     f == :Bx
		real.(getproperty.(x,:A))
	elseif f == :By
		imag.(getproperty.(x,:A))
	elseif f == :Δf
		getproperty.(x,:Δf)
	elseif f == :T || f == :delay
		getproperty.(x[1,:],f)
	elseif f == :dur
		T, delay = x.T, x.delay
		ΔT = [sum(t) for t=T] .+ delay
		ΔT
	else
		getproperty.(x,f)
	end
end

# RF comparison
Base.isapprox(rf1::RF, rf2::RF) = begin
    return all(length(getfield(rf1, k)) == length(getfield(rf2, k)) for k ∈ fieldnames(RF))
        all(≈(getfield(rf1, k), getfield(rf2, k), atol=1e-9) for k ∈ fieldnames(RF))
end

# Properties
size(r::RF, i::Int64) = 1 #To fix [r;r;;] concatenation of Julia 1.7.3
*(α::Complex{T}, x::RF) where {T<:Real} = RF(α*x.A,x.T,x.Δf,x.delay)

"""
    time = dur(rf::RF)
    time = dur(rf::Array{RF,1})
    time = dur(rf::Array{RF,2})

Duration time in seconds of an RF struct or RF array.

# Arguments
- `rf`: (`::RF` or `::Array{RF,1}` or `::Array{RF,2}`) RF struct or RF array

# Returns
- `time`: (`::Real`, [`s`]) duration of the RF struct or RF array
"""
dur(x::RF) = sum(x.T)
dur(x::Array{RF,1}) = sum(sum(x[i].T) for i=1:size(x,1))
dur(x::Array{RF,2}) = maximum(sum([sum(x[i,j].T) for i=1:size(x,1),j=1:size(x,2)],dims=2))

"""
    rf = RF(f::Function, T::Real, N::Int64)

Generates an RF sequence with amplitudes sampled from a function waveform.

# Arguments
- `f`: (`::Function`, [`T`]) function representing the RF amplitude waveform
- `T`: (`::Real`, [`s`]) duration of the RF pulse
- `N`: (`::Int64`) number of samples in the RF pulse

# Returns
- `rf`:(`::RF`) RF struct with amplitude defined by the function `f`
"""
RF(f::Function, T::Real, N::Int64=301; delay::Real=0, Δf=0) = begin
	t = range(0,T;length=N)
	A = f.(t)
	RF(A,T,Δf,delay)
end

"""
    α = get_flip_angle(rf::RF)

Calculates the flip angle α [deg] of an RF struct. α = γ ∫ B1(τ) dτ

# Arguments
- `rf`: (`::RF`) RF struct

# Returns
- `α`: (`::Real`, `[deg]`) flip angle RF struct `rf`
"""
get_flip_angle(x::RF) = begin
	A, NA, T, NT = x.A, length(x.A), x.T, length(x.T)
	dT = T / NA * NT
	α = round(360 * γ * abs(sum(A .* dT)), digits=3) #Pulseq
	return α
end

"""
    time = get_RF_center(rf::RF)

Calculates the time where is the center of the RF pulse `rf`. This calculation includes the
RF delay.

# Arguments
- `rf`: (`::RF`) RF struct

# Returns
- `time`: (`::Real`, `[s]`) time where is the center of the RF pulse `rf`
"""
get_RF_center(x::RF) = begin
	A, NA, T, NT, delay = x.A, length(x.A), x.T, length(x.T), x.delay
	dT = T / NA * NT .* ones(NA)
	t = cumsum([0; dT])[1:end-1]
	t_center = sum(abs.(A) .* t) ./ sum(abs.(A))
	idx = argmin(abs.(t .- t_center))
	t_center += delay + dT[idx]/2
	return t_center
end

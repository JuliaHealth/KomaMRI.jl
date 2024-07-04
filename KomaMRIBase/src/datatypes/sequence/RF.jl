"""
    rf = RF(A, T)
    rf = RF(A, T, Δf)
    rf = RF(A, T, Δf, delay)

The RF struct represents a Radio Frequency excitation of a sequence event.

# Arguments
- `A`: (`::Complex`, `[T]`) RF complex amplitud modulation (AM), ``B_1(t) = |B_1(t)|
    e^{i\\phi(t)} = B_{1}(t) + iB_{1,y}(t) ``
- `T`: (`::Real`, [`s`]) RF duration
- `Δf`: (`::Real` or `::Vector`, [`Hz`]) RF frequency difference with respect to the Larmor frequency.
    This can be a number but also a vector to represent frequency modulated signals (FM).
- `delay`: (`::Real`, [`s`]) RF delay time

# Returns
- `rf`: (`::RF`) the RF struct

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
        return if any(T .< 0) || delay < 0
            error("RF timings must be non-negative.")
        else
            new(A, T, Δf, delay)
        end
    end
    function RF(A, T, Δf)
        return any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, Δf, 0.0)
    end
    function RF(A, T)
        return any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, 0.0, 0.0)
    end
end

"""
    str = show(io::IO, x::RF)

Displays information about the RF struct `x` in the julia REPL.

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `str`: (`::String`) output string message
"""
Base.show(io::IO, x::RF) = begin
    r(x) = round.(x, digits=4)
    compact = get(io, :compact, false)
    if !compact
        wave = length(x.A) == 1 ? r(x.A * 1e6) : "∿"
        print(
            io,
            (x.delay > 0 ? "←$(r(x.delay*1e3)) ms→ " : "") *
            "RF($(wave) uT, $(r(sum(x.T)*1e3)) ms, $(r(x.Δf)) Hz)",
        )
    else
        wave = length(x.A) == 1 ? "⊓" : "∿"
        print(io, (sum(abs.(x.A)) > 0 ? wave : "⇿") * "($(r((x.delay+sum(x.T))*1e3)) ms)")
    end
end

"""
    y = getproperty(x::Vector{RF}, f::Symbol)
    y = getproperty(x::Matrix{RF}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the RF vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{RF}` or `::Matrix{RF}`) vector or matrix of RF structs
- `f`: (`::Symbol`, opts: [`:A`, `:Bx`, `:By`, `:T`, `:Δf`, `:delay` and `:dur`]) input
    symbol that represents a property of the vector or matrix of RF structs

# Returns
- `y`: (`::Vector{Any}` or `::Matrix{Any}`) vector with the property defined by the
    symbol `f` for all elements of the RF vector or matrix `x`
"""
getproperty(x::Vector{RF}, f::Symbol) = getproperty.(x, f)
getproperty(x::Matrix{RF}, f::Symbol) = begin
    if f == :Bx
        real.(getproperty.(x, :A))
    elseif f == :By
        imag.(getproperty.(x, :A))
    elseif f == :Δf
        getproperty.(x, :Δf)
    elseif f == :T || f == :delay
        getproperty.(x[1, :], f)
    elseif f == :dur
        T, delay = x.T, x.delay
        ΔT = [sum(t) for t in T] .+ delay
        ΔT
    else
        getproperty.(x, f)
    end
end

"""
For comparing two `RF`s custom types
"""
function Base.isapprox(rf1::RF, rf2::RF)
    for k in fieldnames(RF)
        if length(getfield(rf1, k)) != length(getfield(rf2, k))
            return false
        end
    end
    return all(≈(getfield(rf1, k), getfield(rf2, k), atol=1e-9) for k in fieldnames(RF))
end

# Properties
size(r::RF, i::Int64) = 1 #To fix [r;r;;] concatenation of Julia 1.7.3
*(α::Complex{T}, x::RF) where {T<:Real} = RF(α * x.A, x.T, x.Δf, x.delay)
*(α::Real, x::RF) = RF(α * x.A, x.T, x.Δf, x.delay)

"""
    y = dur(x::RF)
    y = dur(x::Array{RF,1})
    y = dur(x::Array{RF,2})

Duration time in [s] of RF struct or RF array.

# Arguments
- `x`: (`::RF` or `::Array{RF,1}` or `::Array{RF,2}`) RF struct or RF array

# Returns
- `y`: (`::Float64`, [`s`]) duration of the RF struct or RF array
"""
dur(x::RF) = sum(x.T)
dur(x::Array{RF,1}) = sum(sum(x[i].T) for i in 1:size(x, 1))
function dur(x::Array{RF,2})
    return maximum(sum([sum(x[i, j].T) for i in 1:size(x, 1), j in 1:size(x, 2)]; dims=2))
end

"""
    rf = RF_fun(f::Function, T::Real, N::Int64)

Generate an RF sequence with amplitudes sampled from a function waveform.

!!! note
    This function is not being used in this KomaMRI version.

# Arguments
- `f`: (`::Function`, [`T`]) function for the RF amplitud waveform
- `T`: (`::Real`, [`s`]) duration of the RF pulse
- `N`: (`::Int64`) number of samples of the RF pulse

# Returns
- `rf`:(`::RF`) RF struct with amplitud defined by the function `f`
"""
RF(f::Function, T::Real, N::Int64=301; delay::Real=0, Δf=0) = begin
    t = range(0, T; length=N)
    A = f.(t)
    RF(A, T, Δf, delay)
end

"""
    α = get_flip_angle(x::RF)

Calculates the flip angle α [deg] of an RF struct. α = γ ∫ B1(τ) dτ

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `α`: (`::Int64`, `[deg]`) flip angle RF struct `x`
"""
get_flip_angle(x::RF) = begin
    dt = diff(times(x))
    B1 = ampls(x)
    α = round(360.0 * γ * abs(trapz(dt, B1)); digits=3)
    return α
end

"""
    t = get_RF_center(x::RF)

Calculates the time where is the center of the RF pulse `x`. This calculation includes the
RF delay.

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `t`: (`::Int64`, `[s]`) time where is the center of the RF pulse `x`
"""
get_RF_center(x::RF) = begin
    t = times(x)
    B1 = ampls(x)
    t_center = sum(abs.(B1) .* t) ./ sum(abs.(B1))
    return t_center
end

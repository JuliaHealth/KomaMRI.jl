"""
    Rx = rotx(θ::Real)

Rotates vector counter-clockwise with respect to the x-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rx`: (`::Matrix{Int64}`) rotation matrix
"""
rotx(θ::Real) = [
    1  0  0
    0 cos(θ) -sin(θ)
    0 sin(θ) cos(θ)
]

"""
    Ry = roty(θ::Real)

Rotates vector counter-clockwise with respect to the y-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Ry`: (`::Matrix{Int64}`) rotation matrix
"""
roty(θ::Real) = [
    cos(θ) 0 sin(θ)
    0  1  0
    -sin(θ) 0  cos(θ)
]

"""
    Rz = rotz(θ::Real)

Rotates vector counter-clockwise with respect to the z-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rz`: (`::Matrix{Int64}`) rotation matrix
"""
rotz(θ::Real) = [
    cos(θ) -sin(θ) 0
    sin(θ) cos(θ) 0
    0   0   1
]

"""
    gr = Grad(A, T)
    gr = Grad(A, T, rise)
    gr = Grad(A, T, rise, delay)
    gr = Grad(A, T, rise, fall, delay)
    gr = Grad(A, T, rise, fall, delay, first, last)

The Grad struct represents a gradient of a sequence event.

# Arguments
- `A`: (`::Real` or `::Vector`, `[T/m]`) amplitude of the gradient
- `T`: (`::Real` or `::Vector`, `[s]`) duration of the flat-top
- `rise`: (`::Real`, `[s]`) duration of the rise
- `fall`: (`::Real`, `[s]`) duration of the fall
- `delay`: (`::Real`, `[s]`) duration of the delay

# Returns
- `gr`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> gr = Grad(1, 1, 0.1, 0.1, 0.2)

julia> seq = Sequence([gr]); plot_seq(seq)
```
"""
mutable struct Grad
    A
    T
    rise::Real
    fall::Real
    delay::Real
    first
    last
    function Grad(A, T, rise, fall, delay)
        return if all(T .< 0) || rise < 0 || fall < 0 || delay < 0
            error("Gradient timings must be positive.")
        else
            new(A, T, rise, fall, delay, 0.0, 0.0)
        end
    end
    function Grad(A, T, rise, delay)
        return if all(T .< 0) < 0 || rise < 0 || delay < 0
            error("Gradient timings must be positive.")
        else
            new(A, T, rise, rise, delay, 0.0, 0.0)
        end
    end
    function Grad(A, T, rise)
        return if all(T .< 0) < 0 || rise < 0
            error("Gradient timings must be positive.")
        else
            new(A, T, rise, rise, 0.0, 0.0, 0.0)
        end
    end
    function Grad(A, T)
        return if all(T .< 0) < 0
            error("Gradient timings must be positive.")
        else
            new(A, T, 0.0, 0.0, 0.0, 0.0, 0.0)
        end
    end
end

"""
    gr = Grad(f::Function, T::Real, N::Integer; delay::Real)

Generates an arbitrary gradient waveform defined by the function `f` in the interval t ∈
[0,`T`]. The time separation between two consecutive samples is given by T/(N-1).

# Arguments
- `f`: (`::Function`) function that describes the gradient waveform
- `T`: (`::Real`, `[s]`) duration of the gradient waveform
- `N`: (`::Integer`, `=300`) number of samples of the gradient waveform

# Keywords
- `delay`: (`::Real`, `=0`, `[s]`) delay time of the waveform

# Returns
- `gr`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> gx = Grad(t -> sin(π*t / 0.8), 0.8)

julia> seq = Sequence([gx]); plot_seq(seq)
```
"""
Grad(f::Function, T::Real, N::Integer=300; delay::Real=0) = begin
    t = range(0.0, T; length=N)
    G = f.(t)
    return Grad(G, T, 0.0, 0.0, delay)
end

"""
    str = show(io::IO, x::Grad)

Displays information about the Grad struct `x` in the julia REPL.

# Arguments
- `x`: (`::Grad`) Grad struct

# Returns
- `str` (`::String`) output string message
"""
Base.show(io::IO, x::Grad) = begin
    r(x) = round.(x, digits=4)
    compact = get(io, :compact, false)
    if !compact
        wave = length(x.A) == 1 ? r(x.A * 1e3) : "∿"
        if x.rise == x.fall == 0.0
            print(
                io,
                (x.delay > 0 ? "←$(r(x.delay*1e3)) ms→ " : "") *
                "Grad($(wave) mT, $(r(sum(x.T)*1e3)) ms)",
            )
        else
            print(
                io,
                (x.delay > 0 ? "←$(r(x.delay*1e3)) ms→ " : "") *
                "Grad($(wave) mT, $(r(sum(x.T)*1e3)) ms, ↑$(r(x.rise*1e3)) ms, ↓$(r(x.fall*1e3)) ms)",
            )
        end
    else
        wave = length(x.A) == 1 ? "⊓" : "∿"
        print(
            io,
            (is_on(x) ? wave : "⇿") * "($(r((x.delay+x.rise+x.fall+sum(x.T))*1e3)) ms)",
        )
    end
end

"""
    y = getproperty(x::Vector{Grad}, f::Symbol)
    y = getproperty(x::Matrix{Grad}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the Grad vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{Grad}` or `::Matrix{Grad}`) vector or matrix of Grad structs
- `f`: (`::Symbol`, opts: [`:x`, `:y`, `:z`, `:T`, `:delay`, `:rise`, `:delay`, `:dur`,
    `:A`, `f`]) input symbol that represents a property of the vector or matrix of Grad
    structs

# Returns
- `y`: (`::Vector{Any}` or `::Matrix{Any}`) vector or matrix with the property defined
    by the symbol `f` for all elements of the Grad vector or matrix `x`
"""
getproperty(x::Vector{Grad}, f::Symbol) = getproperty.(x, f)
getproperty(x::Matrix{Grad}, f::Symbol) = begin
    if f == :x
        x[1, :]
    elseif f == :y && size(x, 1) >= 2
        x[2, :]
    elseif f == :z && size(x, 1) >= 3
        x[3, :]
    elseif f == :dur
        maximum(dur.(x); dims=1)[:]
    else
        getproperty.(x, f)
    end
end

"""
For comparing two `Grad`s custom types
"""
function Base.isapprox(gr1::Grad, gr2::Grad)
    for k in fieldnames(Grad)
        if length(getfield(gr1, k)) != length(getfield(gr2, k))
            return false
        end
    end
    return all(getfield(gr1, k) ≈ getfield(gr2, k) for k in fieldnames(Grad))
end

# Gradient operations
*(x::Grad, α::Real) = Grad(α * x.A, x.T, x.rise, x.fall, x.delay)
*(α::Real, x::Grad) = Grad(α * x.A, x.T, x.rise, x.fall, x.delay)
function *(A::Matrix{Float64}, GR::Matrix{Grad})
    N, M = size(GR)
    return [sum(A[i, 1:N] .* GR[:, j]) for i in 1:N, j in 1:M]
end
Base.zero(::Grad) = Grad(0.0, 0.0)
size(g::Grad, i::Int64) = 1 #To fix [g;g;;] concatenation of Julia 1.7.3
/(x::Grad, α::Real) = Grad(x.A / α, x.T, x.rise, x.fall, x.delay)
+(x::Grad, y::Grad) = Grad(x.A .+ y.A, x.T, x.rise, x.fall, x.delay)
+(x::Array{Grad,1}, y::Array{Grad,1}) = [x[i] + y[i] for i in 1:length(x)]
-(x::Grad) = -1 * x
-(x::Grad, y::Grad) = Grad(x.A .- y.A, x.T, x.rise, x.fall, x.delay)

# Gradient functions
function vcat(x::Array{Grad,1}, y::Array{Grad,1})
    return [i == 1 ? x[j] : y[j] for i in 1:2, j in 1:length(x)]
end
function vcat(x::Array{Grad,1}, y::Array{Grad,1}, z::Array{Grad,1})
    return [
        if i == 1
            x[j]
        elseif i == 2
            y[j]
        else
            z[j]
        end for i in 1:3, j in 1:length(x)
    ]
end

"""
    y = dur(x::Grad)
    y = dur(x::Vector{Grad})

Duration time in [s] of Grad struct or Grad array. When the input is a gradient vector, then
the duration is the maximum duration of all the elements of the gradient vector.

# Arguments
- `x`: (`::Grad` or `::Vector{Grad}`) RF struct or RF array

# Returns
- `y`: (`::Float64`, `[s]`) duration of the RF struct or RF array
"""
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
dur(x::Vector{Grad}) = maximum(dur.(x); dims=1)[:]
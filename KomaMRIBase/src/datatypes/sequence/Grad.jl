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
mutable struct Grad{AT,TT}
    A::AT
    T::TT
    rise::Float64
    fall::Float64
    delay::Float64
    first::Float64
    last::Float64
    Grad(A, T, rise, fall, delay, first, last) =
        _has_negative_timings(T) || rise < 0 || fall < 0 || delay < 0 ?
        error("Gradient timings must be positive.") :
        new{typeof(A),typeof(T)}(A, T, rise, fall, delay, first, last)
    Grad(A, T, rise, fall, delay)              = Grad(A, T, rise, fall, delay, 0.0, 0.0)
    Grad(A, T, rise, delay)                    = Grad(A, T, rise, rise, delay, 0.0, 0.0)
    Grad(A, T, rise)                           = Grad(A, T, rise, rise, 0.0,   0.0, 0.0)
    Grad(A, T)                                 = Grad(A, T, 0.0,  0.0,  0.0,   0.0, 0.0)
end

const TrapezoidalGrad = Grad{AT,TT} where {AT<:Number,TT<:Number}
const UniformlySampledGrad = Grad{AT,TT} where {AT<:AbstractVector{<:Number},TT<:Number}
const TimeShapedGrad = Grad{AT,TT} where {AT<:AbstractVector{<:Number},TT<:AbstractVector{<:Number}}

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
    y = getproperty(x::Array{<:Grad}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the Grad vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Array{<:Grad}`) vector or matrix of Grad structs
- `f`: (`::Symbol`, opts: [`:x`, `:y`, `:z`, `:T`, `:delay`, `:rise`, `:delay`, `:dur`,
    `:A`, `f`]) input symbol that represents a property of the vector or matrix of Grad
    structs

# Returns
- `y`: (`::Array{Any}`) vector or matrix with the property defined
    by the symbol `f` for all elements of the Grad vector or matrix `x`
"""
getproperty(x::AbstractVecOrMat{<:Grad}, f::Symbol) = begin
    if f == :x
        @view x[1, :]
    elseif f == :y && size(x, 1) >= 2
        @view x[2, :]
    elseif f == :z && size(x, 1) >= 3
        @view x[3, :]
    elseif f == :dur
        dur(x)
    elseif f in fieldnames(Grad)
        getfield.(x, f)
    else
        getfield(x, f)
    end
end

# Gradient comparison
field_isapprox(gr1::Grad, gr2::Grad; kwargs...) = isapprox(gr1, gr2; kwargs...)
Base.isapprox(gr1::Grad, gr2::Grad; kwargs...) = fields_isapprox(gr1, gr2; kwargs...)
Base.copy(gr::Grad) = Grad(_deepcopy_fields(gr)...)

# Gradient operations
# zeros(Grad, M, N)
Base.zero(::Type{Grad}) = Grad(0.0, 0.0)
# Rotation
Base.zero(::Grad) = Grad(0.0, 0.0)
*(α::Real, x::Grad) = Grad(α * x.A, copy(x.T), x.rise, x.fall, x.delay, α * x.first, α * x.last)

function *(A::AbstractMatrix{<:Real}, gr::AbstractVector{<:Grad})
    size(A, 2) == length(gr) || throw(DimensionMismatch("matrix columns must match gradient axes"))
    out = Vector{Grad}(undef, size(A, 1))
    for row in 1:size(A, 1)
        g = Grad(0.0, 0.0)
        for axis in 1:size(A, 2)
            c = A[row, axis]
            iszero(c) && continue
            g += c * gr[axis]
        end
        out[row] = g
    end
    return _concrete_event_array_copy(out)
end

function *(A::AbstractMatrix{<:Real}, GR::AbstractMatrix{<:Grad})
    size(A, 2) == size(GR, 1) || throw(DimensionMismatch("matrix columns must match gradient axes"))
    out = Matrix{Grad}(undef, size(A, 1), size(GR, 2))
    Threads.@threads for block in 1:size(GR, 2)
        out[:, block] = A * view(GR, :, block)
    end
    return _concrete_event_array_copy(out)
end

_same_grad_timing_fields(x, y) =
    x.T == y.T &&
    x.rise == y.rise &&
    x.fall == y.fall &&
    x.delay == y.delay

_same_grad_timing(x::TrapezoidalGrad, y::TrapezoidalGrad) =
    _same_grad_timing_fields(x, y)
_same_grad_timing(x::UniformlySampledGrad, y::UniformlySampledGrad) =
    length(x.A) == length(y.A) &&
    _same_grad_timing_fields(x, y)
_same_grad_timing(x::TimeShapedGrad, y::TimeShapedGrad) =
    length(x.A) == length(y.A) &&
    _same_grad_timing_fields(x, y)
_same_grad_timing(::Grad, ::Grad) = false

_amplitude_roundoff_tol(A) = 1000 * eps(float(maximum(abs, A)))

function _strictly_increasing_knots!(t)
    for i in 2:length(t)
        t[i] <= t[i - 1] && (t[i] = t[i - 1] + MIN_RISE_TIME)
    end
    return t
end

function _sort_unique_knots!(t)
    sort!(t)
    isempty(t) && return t
    tol = min(MIN_RISE_TIME / 10, 100 * eps(maximum(abs, t)))
    i = 1
    for j in 2:length(t)
        if abs(t[j] - t[i]) > tol
            i += 1
            t[i] = t[j]
        end
    end
    resize!(t, i)
    return t
end

_sample_values_at(samples, t) =
    isempty(samples.t) ? zeros(length(t)) : linear_interpolation(samples..., extrapolation_bc=0.0).(t)

function _simplify_knots(t, A)
    length(t) <= 3 && return t, A
    keep = trues(length(t))
    amplitude_tol = _amplitude_roundoff_tol(A)
    for i in 2:(length(t) - 1)
        Ai = A[i - 1] + (A[i + 1] - A[i - 1]) * ((t[i] - t[i - 1]) / (t[i + 1] - t[i - 1]))
        keep[i] = !isapprox(A[i], Ai; rtol=0, atol=amplitude_tol)
    end
    count(keep) < 3 && (keep[end - 1] = true)
    return t[keep], A[keep]
end

function _grad_from_knots(t, A)
    if isempty(t) || all(iszero, A)
        return Grad(zero(eltype(A)), 0.0)
    end
    t, A = _simplify_knots(t, A)
    if length(t) >= 2
        t[end - 1] = t[end] - t[end - 1] <= MIN_RISE_TIME ? t[end] : t[end - 1] + MIN_RISE_TIME
    end
    length(t) < 3 && return Grad(A[1], t[end] - t[1], 0.0, 0.0, t[1], A[1], A[end])
    if length(t) >= 4 && all(a -> isapprox(a, A[2]; rtol=0, atol=_amplitude_roundoff_tol(A)), A[2:(end - 1)])
        return Grad(A[2], sum(diff(t[2:(end - 1)])), t[2] - t[1], t[end] - t[end - 1], t[1], A[1], A[end])
    end
    return Grad(
        A[2:(end - 1)],
        diff(t[2:(end - 1)]),
        t[2] - t[1],
        t[end] - t[end - 1],
        t[1],
        A[1],
        A[end],
    )
end

function +(x::Grad, y::Grad)
    is_on(x) || return copy(y)
    is_on(y) || return copy(x)
    if _same_grad_timing(x, y)
        return Grad(x.A .+ y.A, copy(x.T), x.rise, x.fall, x.delay, x.first + y.first, x.last + y.last)
    end
    sx = _gradient_interpolation_samples(x)
    sy = _gradient_interpolation_samples(y)
    t = _sort_unique_knots!(vcat(sx.t, sy.t))
    A = _sample_values_at(sx, t) .+ _sample_values_at(sy, t)
    return _grad_from_knots(t, A)
end
# Others
*(x::Grad, α::Real) = α * x
/(x::Grad, α::Real) = Grad(x.A / α, copy(x.T), x.rise, x.fall, x.delay, x.first / α, x.last / α)
-(x::Grad) = -1 * x
-(x::Grad, y::Grad) = x + (-y)

# Gradient functions
function vcat(x::Vector{T}, y::Vector{T}) where {T<:Grad}
    return [i == 1 ? x[j] : y[j] for i in 1:2, j in 1:length(x)]
end
function vcat(x::Vector{T}, y::Vector{T}, z::Vector{T}) where {T<:Grad}
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
    y = dur(x::Vector{<:Grad})
    y = dur(x::Matrix{<:Grad})

Duration time in [s] of Grad struct or Grad Array.

# Arguments
- `x`: (`::Grad` or `::Vector{<:Grad}` or `::Matrix{<:Grad}`) Grad struct or Grad Array

# Returns
- `y`: (`::Float64`, `[s]`) duration of the Grad struct or Grad Array
"""
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
dur(x::AbstractVector{<:Grad}) = maximum(dur.(x); dims=1)[:]
dur(x::AbstractMatrix{<:Grad}) = maximum(dur.(x); dims=1)[:]

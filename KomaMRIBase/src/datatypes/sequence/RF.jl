abstract type RFUse end
struct Excitation <: RFUse end
struct Refocusing <: RFUse end
struct Inversion <: RFUse end
struct Saturation <: RFUse end
struct Preparation <: RFUse end
struct Other <: RFUse end
struct Undefined <: RFUse end

get_RF_use_from_char(::Val{'e'}) = Excitation()
get_RF_use_from_char(::Val{'r'}) = Refocusing()
get_RF_use_from_char(::Val{'i'}) = Inversion()
get_RF_use_from_char(::Val{'s'}) = Saturation()
get_RF_use_from_char(::Val{'p'}) = Preparation()
get_RF_use_from_char(::Val{'o'}) = Other()
get_RF_use_from_char(::Val{'u'}) = Undefined()

get_char_from_RF_use(::Excitation)  = 'e'
get_char_from_RF_use(::Refocusing)  = 'r'
get_char_from_RF_use(::Inversion)   = 'i'
get_char_from_RF_use(::Saturation)  = 's'
get_char_from_RF_use(::Preparation) = 'p'
get_char_from_RF_use(::Other)       = 'o'
get_char_from_RF_use(::Undefined)   = 'u'

"""
    rf = RF(A, T)
    rf = RF(A, T, Δf)
    rf = RF(A, T, Δf, delay)

The RF struct represents a Radio Frequency excitation of a sequence event.

# Arguments
- `A`: (`::Complex`, `[T]`) RF complex amplitud modulation (AM), ``B_1(t) = |B_1(t)|
    e^{i\\phi(t)} = B_{1}(t) + iB_{1,y}(t) ``
- `T`: (`::Real`, `[s]`) RF duration
- `Δf`: (`::Real` or `::Vector`, `[Hz]`) RF frequency difference with respect to the Larmor frequency.
    This can be a number but also a vector to represent frequency modulated signals (FM).
- `delay`: (`::Real`, `[s]`) RF delay time
- `center`: (`::Real`, `[s]`) RF center time
- `ϕ`: (`::Real`, `[rad]`) RF phase at `center`
- `use`: (`::RFUse`) RF use type

# Returns
- `rf`: (`::RF`) the RF struct

# Examples
```julia-repl
julia> rf = RF(1, 1, 0, 0.2)

julia> seq = Sequence(); @addblock seq += rf; plot_seq(seq)
```
"""
mutable struct RF{AT,TT,ΔFT}
    A::AT
    T::TT
    Δf::ΔFT
    delay::Float64
    center::Union{Float64, Nothing}
    ϕ::Float64
    use::RFUse
    RF(A, T, Δf, delay, center, ϕ, use, ::Val{:preserve}) =
        new{typeof(A),typeof(T),typeof(Δf)}(A, T, Δf, delay, center, ϕ, use)
    function RF(A, T, Δf, delay, ::Nothing, ϕ, use)
        if _has_negative_timings(T) || delay < 0
            error("RF timings must be non-negative.")
        end
        rf = RF(A, T, Δf, delay, nothing, mod(ϕ, 2π), use, Val(:preserve))
        return RF(A, T, Δf, delay, get_RF_center(rf), ϕ, use)
    end
    function RF(A, T, Δf, delay, center, ϕ, use)
        if _has_negative_timings(T) || delay < 0
            error("RF timings must be non-negative.")
        end
        rf = RF(A, T, Δf, delay, center, 0.0, use, Val(:preserve))
        Avalues = ampls(rf)
        isempty(Avalues) && return RF(A, T, Δf, delay, center, mod(ϕ, 2π), use, Val(:preserve))
        t = times(rf; separate_closing_knot=false) .- rf.delay
        Interpolations.deduplicate_knots!(t; move_knots=true)
        value = linear_interpolation(t, Avalues, extrapolation_bc=Interpolations.Flat())(center)
        ϕcenter = iszero(abs(value)) ? 0.0 : mod(angle(value), 2π)
        Arel = iszero(ϕcenter) ? A : A .* cis(-ϕcenter)
        return RF(Arel, T, Δf, delay, center, mod(ϕ + ϕcenter, 2π), use, Val(:preserve))
    end
    RF(A, T, Δf, delay, center, use::RFUse) = RF(A, T, Δf, delay; center, use)
    RF(A, T, Δf, delay, center, ϕ)          = RF(A, T, Δf, delay; center, ϕ)
    RF(A, T, Δf, delay, center)             = RF(A, T, Δf, delay; center)
    RF(A, T, Δf)                            = RF(A, T, Δf, 0.0)
    RF(A, T)                                = RF(A, T, 0.0, 0.0)
end

const BlockPulseRF = RF{AT,TT,ΔFT} where {AT<:Number,TT,ΔFT}
const UniformlySampledRF = RF{AT,TT,ΔFT} where {AT<:AbstractVector{<:Number},TT<:Number,ΔFT}
const TimeShapedRF = RF{AT,TT,ΔFT} where {AT<:AbstractVector{<:Number},TT<:AbstractVector{<:Number},ΔFT}
const FrequencyModulatedRF = RF{AT,TT,ΔFT} where {AT,TT,ΔFT<:AbstractVector{<:Number}}
const UniformlySampledFrequencyModulatedRF =
    RF{AT,TT,ΔFT} where {AT,TT<:Number,ΔFT<:AbstractVector{<:Number}}
const TimeShapedFrequencyModulatedRF =
    RF{AT,TT,ΔFT} where {AT,TT<:AbstractVector{<:Number},ΔFT<:AbstractVector{<:Number}}

function RF(A, T, Δf, delay; center=nothing, ϕ=0.0, use=nothing)
    rf = RF(A, T, Δf, delay, center, ϕ, Undefined())
    rf.use = isnothing(use) ? (get_flip_angle(rf) <= 90.01 ? Excitation() : Refocusing()) : use
    return rf
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
        wave = length(x.A) == 1 ? r(abs(x.A) * 1e6) : "∿"
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
    y = getproperty(x::Array{<:RF}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the RF vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Array{<:RF}`) vector or matrix of RF structs
- `f`: (`::Symbol`, opts: [`:A`, `:Bx`, `:By`, `:T`, `:Δf`, `:delay` and `:dur`]) input
    symbol that represents a property of the vector or matrix of RF structs

# Returns
- `y`: (`::Array{Any}`) vector or matrix with the property defined by the
    symbol `f` for all elements of the RF vector or matrix `x`
"""
getproperty(x::AbstractVecOrMat{<:RF}, f::Symbol) = begin
    if f == :Bx
        real.(cis.(getfield.(x, :ϕ)) .* getfield.(x, :A))
    elseif f == :By
        imag.(cis.(getfield.(x, :ϕ)) .* getfield.(x, :A))
    elseif f == :dur
        dur(x)
    elseif f in fieldnames(RF)
        getfield.(x, f)
    else
        getfield(x, f)
    end
end

Base.isapprox(u1::RFUse, u2::RFUse; kwargs...) = u1 == u2

field_isapprox(rf1::RF, rf2::RF; kwargs...) = isapprox(rf1, rf2; kwargs...)
function Base.isapprox(rf1::RF, rf2::RF; kwargs...)
    typeof(rf1) === typeof(rf2) || return false
    return fields_isapprox(rf1, rf2; kwargs...)
end
Base.copy(rf::RF) = RF(_deepcopy_fields(rf)..., Val(:preserve))
    
# Properties
size(r::RF, i::Int64) = 1 #To fix [r;r;;] concatenation of Julia 1.7.3
*(α::Number, x::RF) = is_on(x) ? RF(abs(α) * x.A, copy(x.T), copy(x.Δf), x.delay, x.center, mod(x.ϕ + angle(α), 2π), x.use, Val(:preserve)) : copy(x)
*(x::RF, α::Number) = α * x

"""
    y = dur(x::RF)
    y = dur(x::Vector{<:RF})
    y = dur(x::Matrix{<:RF})

Duration time in [s] of RF struct or RF Array.

# Arguments
- `x`: (`::RF` or `::Vector{<:RF}` or `::Matrix{<:RF}`) RF struct or RF array

# Returns
- `y`: (`::Float64`, [`s`]) duration of the RF struct or RF array
"""
dur(x::RF) = x.delay + sum(x.T)
dur(x::AbstractVector{<:RF}) = maximum(dur.(x); dims=1)[:]
dur(x::AbstractMatrix{<:RF}) = maximum(dur.(x); dims=1)[:]

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

Calculates the time where is the center of the RF pulse `x` .
It does not include the RF delay and uses the weighted average of times by amplitude.

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `t`: (`::Real` or `Nothing`, `[s]`) time where is the center of the RF pulse `x`, or `nothing` if the RF amplitude is zero
"""
function get_RF_center(rf::BlockPulseRF)
    !isnothing(rf.center) && return rf.center
    return sum(rf.T) / 2
end
function get_RF_center(rf::RF)
    !isnothing(rf.center) && return rf.center
    weights = abs.(ampls(rf))
    isempty(weights) && return 0.0
    total = sum(weights)
    iszero(total) && return 0.0
    return sum(weights .* (times(rf; separate_closing_knot=false) .- rf.delay)) / total
end

@doc raw"""
    translation = Translation(time, dx, dy, dz)

Translation motion struct. It produces a linear translation of the phantom.
Its fields are the final displacements in the three axes (dx, dy, dz) 
and the start and end time of the translation.

# Arguments
- `time`: (`::AbstractTimeSpan{T<:Real}`, `[s]`) time scale
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z

# Returns
- `translation`: (`::Translation`) Translation struct

# Examples
```julia-repl
julia> tr = Translation(time=TimeRange(0.0, 0.5), dx=0.01, dy=0.02, dz=0.03)
```
"""
@with_kw struct Translation{T<:Real, TS<:AbstractTimeSpan{T}} <: SimpleMotion{T}
    time      :: TS
    dx         :: T
    dy         :: T
    dz         :: T
end

TranslationX(time::AbstractTimeSpan{T}, dx::T) where {T<:Real} = Translation(time, dx, zero(T), zero(T))
TranslationY(time::AbstractTimeSpan{T}, dy::T) where {T<:Real} = Translation(time, zero(T), dy, zero(T))
TranslationZ(time::AbstractTimeSpan{T}, dz::T) where {T<:Real} = Translation(time, zero(T), zero(T), dz)

function displacement_x!(
    ux::AbstractArray{T},
    motion::Translation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.time)
    ux .= t_unit .* motion.dx
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
    motion::Translation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.time)
    uy .= t_unit .* motion.dy
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    motion::Translation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.time)
    uz .= t_unit .* motion.dz
    return nothing
end


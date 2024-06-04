@doc raw"""
    translation = Translation(t_start, t_end, dx, dy, dz)

Translation motion struct. It produces a linear translation of the phantom.
Its fields are the final displacements in the three axes (dx, dy, dz) 
and the start and end times of the translation.

# Arguments
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z
- `t_start`: (`::Real`, `[s]`) initial time 
- `t_end`: (`::Real`, `[s]`) final time 

# Returns
- `translation`: (`::Translation`) Translation struct

# Examples
```julia-repl
julia> tr = Translation(dx=0.01, dy=0.02, dz=0.03, t_start=0.0, t_end=0.5)
```
"""
@with_kw struct Translation{T<:Real} <: SimpleMotionType{T}
    dx         :: T
    dy         :: T
    dz         :: T
    t_start::T = typeof(dx)(0.0)
    t_end::T   = typeof(dx)(0.0)
    @assert t_end >= t_start "t_end must be major or equal than t_start"
end

function displacement_x(
    motion_type::Translation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    return t_unit .* motion_type.dx
end

function displacement_y(
    motion_type::Translation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    return t_unit .* motion_type.dy
end

function displacement_z(
    motion_type::Translation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    return t_unit .* motion_type.dz
end

times(motion_type::Translation) = begin
    return [motion_type.t_start, motion_type.t_end]
end

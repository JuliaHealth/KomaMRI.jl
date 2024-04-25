@doc raw"""
    translation = Translation(t_start, t_end, dx, dy, dz)

Translation motion struct. It produces a translation of the phantom in the three directions x, y and z.

```math
ux=\left\{\begin{matrix}
0, & t <= t_start\\
\frac{dx}{t_end-t_start}(t-t_start), & t_start < t < t_end\\ 
dx, & t >= t_end
\end{matrix}\right.
```

```math
uy=\left\{\begin{matrix}
0, & t <= t_start\\
\frac{dy}{t_end-t_start}(t-t_start), & t_start < t < t_end\\ 
dy, & t >= t_end
\end{matrix}\right.
```

```math
uz=\left\{\begin{matrix}
0, & t <= t_start\\
\frac{dz}{t_end-t_start}(t-t_start), & t_start < t < t_end\\ 
dz, & t >= t_end
\end{matrix}\right.
```

# Arguments
- `t_start`: (`::Real`, `[s]`) initial time 
- `t_end`: (`::Real`, `[s]`) final time 
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z

# Returns
- `translation`: (`::Translation`) Translation struct


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

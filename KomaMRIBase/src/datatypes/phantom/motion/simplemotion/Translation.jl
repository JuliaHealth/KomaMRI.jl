"""
    translation = Translation(ti, tf, dx, dy, dz)

Translation motion struct. It produces a translation of the phantom in the three directions x, y and z.

```math
ux  = \left\{\begin{matrix}
0, & t <= ti \\
\frac{dx}{tf-ti}(t-ti), & ti < t < tf \\ 
dx, & t >= tf
\end{matrix}\right.
```

```math
uy  = \left\{\begin{matrix}
0, & t <= ti \\
\frac{dy}{tf-ti}(t-ti), & ti < t < tf \\ 
dy, & t >= tf
\end{matrix}\right.
```

```math
uz  = \left\{\begin{matrix}
0, & t <= ti \\
\frac{dz}{tf-ti}(t-ti), & ti < t < tf \\ 
dz, & t >= tf
\end{matrix}\right.
```

# Arguments
- `ti`: (`::Real`, `[s]`) initial time 
- `tf`: (`::Real`, `[s]`) final time 
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z

# Returns
- `translation`: (`::Translation`) Translation struct


"""

@with_kw struct Translation{T<:Real} <: SimpleMotionType{T}
    ti::T = 0.0
    tf::T = 0.0
    dx::T = 0.0
    dy::T = 0.0
    dz::T = 0.0
end

displacement_x(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    return t_unit .* motion_type.dx
end

displacement_y(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    return t_unit .* motion_type.dy
end

displacement_z(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    return t_unit .* motion_type.dz
end

get_range(motion_type::Translation) = begin
    return motion_type.ti, motion_type.tf
end

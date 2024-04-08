@doc raw"""
    rotation = Rotation(ti, tf, pitch, roll, yaw)
 
Rotation motion struct. It produces a rotation of the phantom in the three axes: 
x (pitch), y (roll), and z (yaw)

```math
\begin{align*}
ux  &= cos(\gamma)cos(\beta)x \\ 
    &+ (cos(\gamma)sin(\beta)sin(\alpha) - sin(\gamma)cos(\alpha))y\\ 
    &+ (cos(\gamma)sin(\beta)cos(\alpha) + sin(\gamma)sin(\alpha))z\\
    &- x
\end{align*}
```

```math
\begin{align*}
uy  &= sin(\gamma)cos(\beta)x \\ 
    &+ (sin(\gamma)sin(\beta)sin(\alpha) + cos(\gamma)cos(\alpha))y\\ 
    &+ (sin(\gamma)sin(\beta)cos(\alpha) - cos(\gamma)sin(\alpha))z\\
    &- y
\end{align*}
```

```math
\begin{align*}
uz  &= -sin(\beta)x \\ 
    &+ cos(\beta)sin(\alpha)y\\ 
    &+ cos(\beta)cos(\alpha)z\\
    &- z
\end{align*}
```

where:
    
```math
\alpha  = \left\{\begin{matrix}
0, & t <= ti \\
\frac{pitch}{tf-ti}(t-ti), & ti < t < tf \\ 
pitch, & t >= tf
\end{matrix}\right.
,\qquad
\beta  = \left\{\begin{matrix}
0, & t <= ti \\
\frac{roll}{tf-ti}(t-ti), & ti < t < tf \\ 
roll, & t >= tf
\end{matrix}\right.
,\qquad
\gamma  = \left\{\begin{matrix}
0, & t <= ti \\
\frac{yaw}{tf-ti}(t-ti), & ti < t < tf \\ 
yaw, & t >= tf
\end{matrix}\right.
```

# Arguments
- `ti`: (`::Real`, `[s]`) initial time 
- `tf`: (`::Real`, `[s]`) final time 
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z

# Returns
- `rotation`: (`::Rotation`) Rotation struct

"""

@with_kw struct Rotation{T<:Real} <: SimpleMotionType{T} 
    ti::T
    tf::T
    pitch::T = typeof(ti)(0.0)
    roll::T  = typeof(ti)(0.0)
    yaw::T   = typeof(ti)(0.0)
end 

is_composable(motion_type::Rotation) = true

position_x(motion_type::Rotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    α = t_unit .* (motion_type.pitch * π/180)
    β = t_unit .* (motion_type.roll * π/180)
    γ = t_unit .* (motion_type.yaw * π/180)
    return cos.(γ) .* cos.(β) .* x   +   (cos.(γ) .* sin.(β) .* sin.(α) .- sin.(γ) .* cos.(α)) .* y   +   (cos.(γ) .* sin.(β) .* cos.(α) .+ sin.(γ) .* sin.(α)) .* z
end

position_y(motion_type::Rotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    α = t_unit .* (motion_type.pitch * π/180)
    β = t_unit .* (motion_type.roll * π/180)
    γ = t_unit .* (motion_type.yaw * π/180)
    return sin.(γ) .* cos.(β) .* x   +   (sin.(γ) .* sin.(β) .* sin.(α) .+ cos.(γ) .* cos.(α)) .* y   +   (sin.(γ) .* sin.(β) .* cos.(α) .- cos.(γ) .* sin.(α)) .* z
end

position_z(motion_type::Rotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    α = t_unit .* (motion_type.pitch * π/180)
    β = t_unit .* (motion_type.roll * π/180)
    γ = t_unit .* (motion_type.yaw * π/180)
    return -sin.(β) .* x   +   cos.(β) .* sin.(α) .* y
end

get_time_nodes(motion_type::Rotation) = begin
    return [motion_type.ti, motion_type.tf]
end
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

displacement_x(motion_type::Rotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    α = t_unit .* (motion_type.pitch)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.yaw)
    return cosd.(γ) .* cosd.(β) .* x   +   (cosd.(γ) .* sind.(β) .* sind.(α) .- sind.(γ) .* cosd.(α)) .* y   +   (cosd.(γ) .* sind.(β) .* cosd.(α) .+ sind.(γ) .* sind.(α)) .* z .- x
end

displacement_y(motion_type::Rotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    α = t_unit .* (motion_type.pitch)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.yaw)
    return sind.(γ) .* cosd.(β) .* x   +   (sind.(γ) .* sind.(β) .* sind.(α) .+ cosd.(γ) .* cosd.(α)) .* y   +   (sind.(γ) .* sind.(β) .* cosd.(α) .- cosd.(γ) .* sind.(α)) .* z .- y
end

displacement_z(motion_type::Rotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    α = t_unit .* (motion_type.pitch)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.yaw)
    return -sind.(β) .* x   +   cosd.(β) .* sind.(α) .* y .- z
end

time_nodes(motion_type::Rotation) = begin
    return [motion_type.ti, motion_type.tf]
end
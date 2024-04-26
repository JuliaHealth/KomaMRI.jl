@doc raw"""
    rotation = Rotation(t_start, t_end, pitch, roll, yaw)
 
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
0, & t <= t_start \\
\frac{pitch}{t_end-t_start}(t-t_start), & t_start < t < t_end \\ 
pitch, & t >= t_end
\end{matrix}\right.
,\qquad
\beta  = \left\{\begin{matrix}
0, & t <= t_start \\
\frac{roll}{t_end-t_start}(t-t_start), & t_start < t < t_end \\ 
roll, & t >= t_end
\end{matrix}\right.
,\qquad
\gamma  = \left\{\begin{matrix}
0, & t <= t_start \\
\frac{yaw}{t_end-t_start}(t-t_start), & t_start < t < t_end \\ 
yaw, & t >= t_end
\end{matrix}\right.
```

# Arguments
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z
- `t_start`: (`::Real`, `[s]`) initial time 
- `t_end`: (`::Real`, `[s]`) final time 

# Returns
- `rotation`: (`::Rotation`) Rotation struct

"""
@with_kw struct Rotation{T<:Real} <: SimpleMotionType{T}
    pitch      :: T
    roll       :: T
    yaw        :: T
    t_start::T = typeof(pitch)(0.0)
    t_end      = typeof(pitch)(0.0)
    @assert t_end >= t_start "t_end must be major or equal than t_start"
end

is_composable(motion_type::Rotation) = true

function displacement_x(
    motion_type::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    α = t_unit .* (motion_type.pitch)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.yaw)
    return cosd.(γ) .* cosd.(β) .* x +
           (cosd.(γ) .* sind.(β) .* sind.(α) .- sind.(γ) .* cosd.(α)) .* y +
           (cosd.(γ) .* sind.(β) .* cosd.(α) .+ sind.(γ) .* sind.(α)) .* z .- x
end

function displacement_y(
    motion_type::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return sind.(γ) .* cosd.(β) .* x +
           (sind.(γ) .* sind.(β) .* sind.(α) .+ cosd.(γ) .* cosd.(α)) .* y +
           (sind.(γ) .* sind.(β) .* cosd.(α) .- cosd.(γ) .* sind.(α)) .* z .- y
end

function displacement_z(
    motion_type::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return -sind.(β) .* x + cosd.(β) .* sind.(α) .* y .- z
end

times(motion_type::Rotation) = begin
    return [motion_type.t_start, motion_type.t_end]
end
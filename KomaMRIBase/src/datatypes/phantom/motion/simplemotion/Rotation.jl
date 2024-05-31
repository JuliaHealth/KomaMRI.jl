@doc raw"""
    rotation = Rotation(t_start, t_end, pitch, roll, yaw)
 
Rotation motion struct. It produces a rotation of the phantom in the three axes: 
x (pitch), y (roll), and z (yaw).
We follow the RAS (Right-Anterior-Superior) orientation, 
and the rotations are applied following the right-hand rule (counter-clockwise):

| ![Head Rotation Axis](../assets/head-rotation-axis.svg) |
| ------------------------------------------------------- |

The applied rotation matrix is obtained as follows: 
```math
\begin{equation}
\begin{aligned}
R &= R_z(\alpha) R_y(\beta) R_x(\gamma) \\
  &= \begin{bmatrix}
\cos \alpha & -\sin \alpha & 0 \\
\sin \alpha & \cos \alpha & 0 \\
0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
\cos \beta & 0 & \sin \beta \\
0 & 1 & 0 \\
-\sin \beta & 0 & \cos \beta
\end{bmatrix}
\begin{bmatrix}
1 & 0 & 0 \\
0 & \cos \gamma & -\sin \gamma \\
0 & \sin \gamma & \cos \gamma
\end{bmatrix} \\
  &= \begin{bmatrix}
\cos \alpha \cos \beta & \cos \alpha \sin \beta \sin \gamma - \sin \alpha \cos \gamma & \cos \alpha \sin \beta \cos \gamma + \sin \alpha \sin \gamma \\
\sin \alpha \cos \beta & \sin \alpha \sin \beta \sin \gamma + \cos \alpha \cos \gamma & \sin \alpha \sin \beta \cos \gamma - \cos \alpha \sin \gamma \\
-\sin \beta & \cos \beta \sin \gamma & \cos \beta \cos \gamma
\end{bmatrix}
\end{aligned}
\end{equation}
```

# Arguments
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z
- `t_start`: (`::Real`, `[s]`) initial time 
- `t_end`: (`::Real`, `[s]`) final time 

# Returns
- `rotation`: (`::Rotation`) Rotation struct

# Examples
```julia-repl
julia> rt = Rotation(pitch=15.0, roll=0.0, yaw=20.0, t_start=0.1, t_end=0.5)
```
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
    α = t_unit .* (motion_type.yaw)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.pitch)
    return cosd.(α) .* cosd.(β) .* x +
           (cosd.(α) .* sind.(β) .* sind.(γ) .- sind.(α) .* cosd.(γ)) .* y +
           (cosd.(α) .* sind.(β) .* cosd.(γ) .+ sind.(α) .* sind.(γ)) .* z .- x
end

function displacement_y(
    motion_type::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    α = t_unit .* (motion_type.yaw)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.pitch)
    return sind.(α) .* cosd.(β) .* x +
           (sind.(α) .* sind.(β) .* sind.(γ) .+ cosd.(α) .* cosd.(γ)) .* y +
           (sind.(α) .* sind.(β) .* cosd.(γ) .- cosd.(α) .* sind.(γ)) .* z .- y
end

function displacement_z(
    motion_type::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    α = t_unit .* (motion_type.yaw)
    β = t_unit .* (motion_type.roll)
    γ = t_unit .* (motion_type.pitch)
    return -sind.(β) .* x + 
            cosd.(β) .* sind.(γ) .* y +
            cosd.(β) .* cosd.(γ) .* z .- z
end

times(motion_type::Rotation) = begin
    return [motion_type.t_start, motion_type.t_end]
end
@doc raw"""
    rotation = Rotation(time, pitch, roll, yaw)
 
Rotation motion struct. It produces a rotation of the phantom in the three axes: 
x (pitch), y (roll), and z (yaw).
We follow the RAS (Right-Anterior-Superior) orientation, 
and the rotations are applied following the right-hand rule (counter-clockwise):

![Head Rotation Axis](../assets/head-rotation-axis.svg)

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
- `time`: (`::AbstractTimeSpan{T<:Real}`, `[s]`) time scale
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z

# Returns
- `rotation`: (`::Rotation`) Rotation struct

# Examples
```julia-repl
julia> rt = Rotation(time=TimeRange(0.0, 0.5), pitch=15.0, roll=0.0, yaw=20.0)
```
"""
@with_kw struct Rotation{T<:Real, TS<:AbstractTimeSpan{T}} <: SimpleMotion{T}
    time      :: TS
    pitch      :: T
    roll       :: T
    yaw        :: T
end

RotationX(time::AbstractTimeSpan{T}, pitch::T) where {T<:Real} = Rotation(time, pitch, zero(T), zero(T))
RotationY(time::AbstractTimeSpan{T}, roll::T) where {T<:Real}  = Rotation(time, zero(T), roll, zero(T))
RotationZ(time::AbstractTimeSpan{T}, yaw::T) where {T<:Real}   = Rotation(time, zero(T), zero(T), yaw)

is_composable(motion::Rotation) = true

function displacement_x!(
    ux::AbstractArray{T},
    motion::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.time)
    α = t_unit .* (motion.yaw)
    β = t_unit .* (motion.roll)
    γ = t_unit .* (motion.pitch)
    ux .= cosd.(α) .* cosd.(β) .* x +
         (cosd.(α) .* sind.(β) .* sind.(γ) .- sind.(α) .* cosd.(γ)) .* y +
         (cosd.(α) .* sind.(β) .* cosd.(γ) .+ sind.(α) .* sind.(γ)) .* z .- x
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
    motion::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.time)
    α = t_unit .* (motion.yaw)
    β = t_unit .* (motion.roll)
    γ = t_unit .* (motion.pitch)
    uy .= sind.(α) .* cosd.(β) .* x +
         (sind.(α) .* sind.(β) .* sind.(γ) .+ cosd.(α) .* cosd.(γ)) .* y +
         (sind.(α) .* sind.(β) .* cosd.(γ) .- cosd.(α) .* sind.(γ)) .* z .- y
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    motion::Rotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.time)
    α = t_unit .* (motion.yaw)
    β = t_unit .* (motion.roll)
    γ = t_unit .* (motion.pitch)
    uz .=  -sind.(β) .* x + 
            cosd.(β) .* sind.(γ) .* y +
            cosd.(β) .* cosd.(γ) .* z .- z
    return nothing
end
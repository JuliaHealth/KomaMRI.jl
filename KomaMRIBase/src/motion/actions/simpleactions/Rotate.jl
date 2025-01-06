@doc raw"""
    rotate = Rotate(pitch, roll, yaw)
 
Rotate struct. It produces a rotation in the three axes: 
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
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z

# Returns
- `rotate`: (`::Rotate`) Rotate struct

# Examples
```julia-repl
julia> rotate = Rotate(pitch=15.0, roll=0.0, yaw=20.0)
```
"""
@with_kw struct Rotate{T<:Real} <: SimpleAction{T}
    pitch      :: T
    roll       :: T
    yaw        :: T
end

RotateX(pitch::T) where {T<:Real} = Rotate(pitch, zero(T), zero(T))
RotateY(roll::T) where {T<:Real}  = Rotate(zero(T), roll, zero(T))
RotateZ(yaw::T) where {T<:Real}   = Rotate(zero(T), zero(T), yaw)

is_composable(action::Rotate) = true

function displacement_x!(ux, action::Rotate, x, y, z, t)   
    # Not using sind and cosd functions until bug with oneAPI is solved: 
    # https://github.com/JuliaGPU/oneAPI.jl/issues/65
    α = t .* (action.yaw*π/180)
    β = t .* (action.roll*π/180)
    γ = t .* (action.pitch*π/180)
    ux .= cos.(α) .* cos.(β) .* x +
         (cos.(α) .* sin.(β) .* sin.(γ) .- sin.(α) .* cos.(γ)) .* y +
         (cos.(α) .* sin.(β) .* cos.(γ) .+ sin.(α) .* sin.(γ)) .* z .- x
    return nothing
end

function displacement_y!(uy, action::Rotate, x, y, z, t)
    α = t .* (action.yaw*π/180)
    β = t .* (action.roll*π/180)
    γ = t .* (action.pitch*π/180)
    uy .= sin.(α) .* cos.(β) .* x +
         (sin.(α) .* sin.(β) .* sin.(γ) .+ cos.(α) .* cos.(γ)) .* y +
         (sin.(α) .* sin.(β) .* cos.(γ) .- cos.(α) .* sin.(γ)) .* z .- y
    return nothing
end

function displacement_z!(uz, action::Rotate, x, y, z, t)
    α = t .* (action.yaw*π/180)
    β = t .* (action.roll*π/180)
    γ = t .* (action.pitch*π/180)
    uz .=  -sin.(β) .* x + 
            cos.(β) .* sin.(γ) .* y +
            cos.(β) .* cos.(γ) .* z .- z
    return nothing
end
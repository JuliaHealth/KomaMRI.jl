struct CenterOfMass end
Base.:(≈)(::CenterOfMass, ::CenterOfMass) = true
Base.:(≈)(a::CenterOfMass, b) = false
Base.:(≈)(a, b::CenterOfMass) = false

@doc raw"""
    r = Rotate(pitch, roll, yaw, center=CenterOfMass())
 
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
- `center`: (`::NTuple{3,Real}` or `::CenterOfMass`) optional center of rotation, given in global coordinates. Default is center of mass.

# Notes
- Rotations are applied around the point specified in `center`. If omitted, the rotation is centered at the phantom’s center of mass.
- If `center` is not `::CenterOfMass`, the rotation center is interpreted as a fixed point in space (absolute/global coordinates).
- This design ensures that consecutive or inverse rotations behave consistently and predictably, since the rotation center does not change with object transformations.

# Returns
- `r`: (`::Rotate`) Rotate struct

# Examples
```julia-repl
julia> r = Rotate(pitch=15.0, roll=0.0, yaw=20.0)

julia> r = Rotate(pitch=0.0, roll=45.0, yaw=0.0, center=(5e-3,0.0,0.0))
# Rotates around a point 5 mm to the right of the center of mass
```
"""
@with_kw struct Rotate{T<:Real} <: SimpleAction{T}
    pitch      :: T
    roll       :: T
    yaw        :: T
    center     :: Union{CenterOfMass,NTuple{3,T}} = CenterOfMass()
end

RotateX(pitch::T) where {T<:Real} = Rotate(pitch=pitch,   roll=zero(T), yaw=zero(T))
RotateY(roll::T)  where {T<:Real} = Rotate(pitch=zero(T), roll=roll,    yaw=zero(T))
RotateZ(yaw::T)   where {T<:Real} = Rotate(pitch=zero(T), roll=zero(T), yaw=yaw)

get_center(center::CenterOfMass, x, y, z) = (sum(x) / length(x), sum(y) / length(y), sum(z) / length(z))
get_center(center::NTuple, x, y, z)       = center

function displacement_x!(ux, action::Rotate, x, y, z, t)
    cx, cy, cz = get_center(action.center, x, y, z)
    α = deg2rad(action.yaw)
    β = deg2rad(action.roll)
    γ = deg2rad(action.pitch)
    
    @. ux =  cos(α * t) * cos(β * t) * (x - cx) +
            (cos(α * t) * sin(β * t) * sin(γ * t) - sin(α * t) * cos(γ * t)) * (y - cy) +
            (cos(α * t) * sin(β * t) * cos(γ * t) + sin(α * t) * sin(γ * t)) * (z - cz) + 
            cx - x
            
    return nothing
end

function displacement_y!(uy, action::Rotate, x, y, z, t)
    cx, cy, cz = get_center(action.center, x, y, z)
    α = deg2rad(action.yaw)
    β = deg2rad(action.roll)
    γ = deg2rad(action.pitch)
    
    @. uy =  sin(α * t) * cos(β * t) * (x - cx) +
            (sin(α * t) * sin(β * t) * sin(γ * t) + cos(α * t) * cos(γ * t)) * (y - cy) +
            (sin(α * t) * sin(β * t) * cos(γ * t) - cos(α * t) * sin(γ * t)) * (z - cz) + 
            cy - y
            
    return nothing
end

function displacement_z!(uz, action::Rotate, x, y, z, t)
    cx, cy, cz = get_center(action.center, x, y, z)
    α = deg2rad(action.yaw)
    β = deg2rad(action.roll)
    γ = deg2rad(action.pitch)
    
    @. uz = -sin(β * t) * (x - cx) + 
             cos(β * t) * sin(γ * t) * (y - cy) +
             cos(β * t) * cos(γ * t) * (z - cz) + 
             cz - z
            
    return nothing
end
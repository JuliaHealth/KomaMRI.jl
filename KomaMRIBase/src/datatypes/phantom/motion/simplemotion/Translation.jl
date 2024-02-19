"""
Simple Translation Movement - Constant Velocity

Parameters:
- Direction
- Velocity
"""

mutable struct Translation{T<:Real} <: SimpleMotionType{T}
    direction::AbstractVector{T}     # Translation direction 
    v::T                             # Velocity [m/s]
end

function SimpleMotion(type::Translation)
    # Normalize direction
    n = normalize(type.direction)

    # Obtain velocity vector
    v = type.v*n

    ux(x,y,z,t) = v[1].*t
    uy(x,y,z,t) = v[2].*t
    uz(x,y,z,t) = v[3].*t  

    return SimpleMotion(type,ux,uy,uz)
end
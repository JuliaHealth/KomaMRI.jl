"""
Simple Translation Movement - Constant Velocity

Parameters:
- Direction
- Velocity
"""

@with_kw mutable struct Translation{T<:Real} <: SimpleMotionType
    direction::AbstractVector{T}     # Translation direction 
    v::T                             # Velocity [m/s]
end

function SimpleMotion(trans::Translation)
    # Normalize direction
    n = normalize(trans.direction)

    # Obtain velocity vector
    v = trans.v*n

    ux(x,y,z,t) = @view(v[1]).*t
    uy(x,y,z,t) = @view(v[2]).*t
    uz(x,y,z,t) = @view(v[3]).*t  

    return SimpleMotion(type=trans,ux=ux,uy=uy,uz=uz)
end
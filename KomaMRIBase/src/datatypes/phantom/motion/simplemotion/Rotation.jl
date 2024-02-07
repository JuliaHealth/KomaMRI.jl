"""
Simple Rotation Movement

Parameters:
- Rotation axis
- Angular velocity
"""

@with_kw mutable struct Rotation{T<:Real} <: SimpleMotionType
    axis::AbstractVector{T}     # Rotation axis vector
    point::AbstractVector{T}    # Rotation axis point
    f::T                        # Angular velocity [Hz]
end

function get_simple_motion(rot::Rotation)
    # Calculate rotation matrix
    R(t) =  axis_angle(rot.axis, 2π*rot.f*t)
    displace(x,y,z) =  [[x,y,z][i] .- rot.point[i] for i in 1:length(rot.point)]

    ux(x,y,z,t) = begin
        xd, yd, zd = displace(x,y,z)
        xr = hcat([(R(ti)[1,1]*xd + R(ti)[1,2]*yd + R(ti)[1,3]*zd) .+ rot.point[1] - x for ti in t]...)
        return xr
    end

    return SimpleMotion(ux=ux)
end
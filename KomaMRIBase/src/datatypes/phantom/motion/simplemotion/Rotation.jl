"""
Simple Rotation Movement

Parameters:
- Rotation axis
- Angular velocity
"""

struct Rotation{T<:Real} <: SimpleMotionType
    axis::AbstractVector{T}     # Rotation axis vector
    point::AbstractVector{T}    # Rotation axis point
    ω::T                        # Angular velocity [Hz]
end

function get_simple_motion(type::Rotation)
    # return SimpleMotion(
    #     ux = sin(type.v * t)
    # )
end
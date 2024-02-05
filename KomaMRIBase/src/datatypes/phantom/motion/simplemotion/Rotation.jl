"""
Simple Rotation Movement

Parameters:
- Rotation axis
- Angular velocity
"""

struct Rotation <: SimpleMotionType
# 
end

function get_simple_motion(type::Rotation)
    # return SimpleMotion(
    #     ux = sin(type.v * t)
    # )
end
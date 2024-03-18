"""
No Motion

x = x
"""

struct NoMotion <: MotionModel end

Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                q::Union{AbstractRange,AbstractVector,Colon}) = motion
function get_spin_coords(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where T<:Real
    return x, y, z
end

function get_range(motion::NoMotion)
    return 0.0
end
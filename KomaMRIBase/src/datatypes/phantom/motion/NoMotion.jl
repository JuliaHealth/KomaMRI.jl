"""
No Motion

x = x
"""

struct NoMotion{T<:Real} <: MotionModel{T} end

Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}, q::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.:(==)(m1::NoMotion, m2::NoMotion) = true

function get_spin_coords(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where T<:Real
    return x, y, z
end

function get_times(motion::NoMotion)
    return [0.0]
end
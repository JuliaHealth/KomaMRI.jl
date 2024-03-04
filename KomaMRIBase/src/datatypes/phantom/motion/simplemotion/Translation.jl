"""
Simple Translation Movement - Constant Velocity

Parameters:
- Offset
- Velocity

ux = x0 + vt
"""

@with_kw mutable struct Translation{T<:Real} <: SimpleMotionType{T}
    offset::Vector{T}   = zeros(3)   # offset
    velocity::Vector{T} = zeros(3)   # Velocity [m/s]
end

displacement_x(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = motion_type.offset[1] .+ motion_type.velocity[1] .* t
displacement_y(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = motion_type.offset[2] .+ motion_type.velocity[2] .* t
displacement_z(motion_type::Translation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = motion_type.offset[3] .+ motion_type.velocity[3] .* t

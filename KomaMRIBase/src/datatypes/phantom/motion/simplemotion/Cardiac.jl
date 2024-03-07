"""
Cardiac Motion
 
Parameters:
- RR            [Nx1], where N is the number of thr different RR interval durations

ux = 0          (Pending)
"""

@with_kw struct Cardiac{T<:Real} <: SimpleMotionType{T}
    RR::Vector{T} = [1.0]
end

displacement_x(motion_type::Cardiac{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = 0
displacement_y(motion_type::Cardiac{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = 0
displacement_z(motion_type::Cardiac{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = 0
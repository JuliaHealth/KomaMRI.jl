"""
Simple Rotation Movement

Parameters:
- Rotation axis
- Angular velocity
"""

@with_kw mutable struct Rotation{T<:Real} <: SimpleMotionType{T}
    axis::AbstractVector{T}  = [.0, .0, 1]      # Rotation axis vector
    point::AbstractVector{T} = zeros(3)         # Rotation axis point
    f::T = 1.0                                  # Angular velocity [Hz]
end

# Rotation matrix
R(axis, f, t) = Un.(2π*f*t, Ref(axis))

# Displace using the rotation point
displace(x, y, z, point) =  [[x,y,z][i] .- point[i] for i in 1:length(point)]

# Conversion of x,y,z into a (3 x Ns) matrix
matrix(x,y,z) = hcat(x,y,z)'

rotation(x, y, z, t, motion_type::Rotation) = begin
    xd, yd, zd = displace(x, y, z, motion_type.point)
    xyz = matrix(xd, yd, zd)
    rotated = map(x -> x*xyz, R(motion_type.axis, motion_type.f, t))
    return cat(rotated...,dims=3)
end

displacement_x(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = rotation(x,y,z,t,motion_type)[1,:,:] .+ (motion_type.point[1] .- x) 
displacement_y(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = rotation(x,y,z,t,motion_type)[2,:,:] .+ (motion_type.point[2] .- y) 
displacement_z(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = rotation(x,y,z,t,motion_type)[3,:,:] .+ (motion_type.point[3] .- z) 

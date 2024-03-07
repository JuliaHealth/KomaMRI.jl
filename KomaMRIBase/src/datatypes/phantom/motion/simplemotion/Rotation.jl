"""
Rotation Motion
 
Parameters:
- Offset            (α) [1x1] 
- Rotation axis     (n) [3x1]
- Rotation center   (c) [3x1]
- Angular velocity  (ω) [1x1]

ux = Un(ωt)(x - c) + (c - x)

where Un(wt) = Icos(ωt) + sin(ωt)cross(n) + (1-cos(ωt))n.n'
"""

@with_kw struct Rotation{T<:Real} <: SimpleMotionType{T}
    offset::T = 0.0                          # Initial rotation [rad]
    rotation_axis::Vector{T}  = [0, 0, 1.0]  # Rotation axis vector
    rotation_center::Vector{T} = zeros(3)    # Rotation axis point  (TODO: check if this works correctly)
    angular_velocity::T = 2*π                # Angular velocity [rad/s]
end

displacement_x(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    x .-= motion_type.rotation_center[1]
    θ = motion_type.offset .+ motion_type.angular_velocity*t 
    nx, ny, nz = motion_type.rotation_axis
    # Rodrigues' formula
    xr = cos.(θ) .* x + sin.(θ) .* (-nz*y + ny*z) + (1 .- cos.(θ)) .* (nx*(nx*x + ny*y + nz*z))
    xr .+= (motion_type.rotation_center[1] .- x)
    return xr
end

displacement_y(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    y .-= motion_type.rotation_center[2]
    θ = motion_type.offset .+ motion_type.angular_velocity*t 
    nx, ny, nz = motion_type.rotation_axis
    # Rodrigues' formula
    yr = cos.(θ) .* y + sin.(θ) .* (nz*x - nx*z)  + (1 .- cos.(θ)) .* (ny*(nx*x + ny*y + nz*z))
    yr .+= (motion_type.rotation_center[2] .- y)
    return yr
end

displacement_z(motion_type::Rotation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    z .-= motion_type.rotation_center[3]
    θ = motion_type.offset .+ motion_type.angular_velocity*t 
    nx, ny, nz = motion_type.rotation_axis
    # Rodrigues' formula
    zr = cos.(θ) .* z + sin.(θ) .* (-ny*x + nx*y) + (1 .- cos.(θ)) .* (nz*(nx*x + ny*y + nz*z))
    zr .+= (motion_type.rotation_center[3] .- z)
    return zr
end
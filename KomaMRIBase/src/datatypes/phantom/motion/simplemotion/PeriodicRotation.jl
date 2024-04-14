@doc raw"""
    periodic_rotation = PeriodicRotation(period, asymmetry, pitch, roll, yaw)
 
PeriodicRotation motion struct. It produces a rotation of the phantom in the three axes: 
x (pitch), y (roll), and z (yaw)

# Arguments
- `period`: (`::Real`, `[s]`) period 
- `asymmetry`: (`::Real`)  asymmetry factor, between 0 and 1
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z

# Returns
- `periodic_rotation`: (`::PeriodicRotation`) Rotation struct

"""

@with_kw struct PeriodicRotation{T<:Real} <: SimpleMotionType{T}
    period::T
    asymmetry::T = typeof(period)(0.5)
    pitch::T     = typeof(period)(0.0)
    roll::T      = typeof(period)(0.0)
    yaw::T       = typeof(period)(0.0)
end 

is_composable(motion_type::PeriodicRotation) = true

displacement_x(motion_type::PeriodicRotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return cosd.(γ) .* cosd.(β) .* x   +   (cosd.(γ) .* sind.(β) .* sind.(α) .- sind.(γ) .* cosd.(α)) .* y   +   (cosd.(γ) .* sind.(β) .* cosd.(α) .+ sind.(γ) .* sind.(α)) .* z .- x
end

displacement_y(motion_type::PeriodicRotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return sind.(γ) .* cosd.(β) .* x   +   (sind.(γ) .* sind.(β) .* sind.(α) .+ cosd.(γ) .* cosd.(α)) .* y   +   (sind.(γ) .* sind.(β) .* cosd.(α) .- cosd.(γ) .* sind.(α)) .* z .- y
end

displacement_z(motion_type::PeriodicRotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return -sind.(β) .* x   +   cosd.(β) .* sind.(α) .* y .- z
end

time_nodes(motion_type::PeriodicRotation) = begin
    return [0, motion_type.period * motion_type.asymmetry, motion_type.period]
end
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

position_x(motion_type::PeriodicRotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* (motion_type.pitch * π/180)
    β = t_unit .* (motion_type.roll * π/180)
    γ = t_unit .* (motion_type.yaw * π/180)
    return cos.(γ) .* cos.(β) .* x   +   (cos.(γ) .* sin.(β) .* sin.(α) .- sin.(γ) .* cos.(α)) .* y   +   (cos.(γ) .* sin.(β) .* cos.(α) .+ sin.(γ) .* sin.(α)) .* z
end

position_y(motion_type::PeriodicRotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* (motion_type.pitch * π/180)
    β = t_unit .* (motion_type.roll * π/180)
    γ = t_unit .* (motion_type.yaw * π/180)
    return sin.(γ) .* cos.(β) .* x   +   (sin.(γ) .* sin.(β) .* sin.(α) .+ cos.(γ) .* cos.(α)) .* y   +   (sin.(γ) .* sin.(β) .* cos.(α) .- cos.(γ) .* sin.(α)) .* z
end

position_z(motion_type::PeriodicRotation{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* (motion_type.pitch * π/180)
    β = t_unit .* (motion_type.roll * π/180)
    γ = t_unit .* (motion_type.yaw * π/180)
    return -sin.(β) .* x   +   cos.(β) .* sin.(α) .* y
end

get_time_nodes(motion_type::PeriodicRotation) = begin
    return [0, motion_type.period * motion_type.asymmetry, motion_type.period]
end
@doc raw"""
    periodic_rotation = PeriodicRotation(period, asymmetry, pitch, roll, yaw)
 
PeriodicRotation motion struct. It produces a rotation of the phantom in the three axes: 
x (pitch), y (roll), and z (yaw)

# Arguments
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z
- `period`: (`::Real`, `[s]`) period 
- `asymmetry`: (`::Real`)  asymmetry factor, between 0 and 1

# Returns
- `periodic_rotation`: (`::PeriodicRotation`) PeriodicRotation struct

"""
@with_kw struct PeriodicRotation{T<:Real} <: SimpleMotionType{T}
    pitch        :: T
    roll         :: T
    yaw          :: T
    period::T    = typeof(pitch)(0.0)
    asymmetry::T = typeof(pitch)(0.5)
end

is_composable(motion_type::PeriodicRotation) = true

function displacement_x(
    motion_type::PeriodicRotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return cosd.(γ) .* cosd.(β) .* x +
           (cosd.(γ) .* sind.(β) .* sind.(α) .- sind.(γ) .* cosd.(α)) .* y +
           (cosd.(γ) .* sind.(β) .* cosd.(α) .+ sind.(γ) .* sind.(α)) .* z .- x
end

function displacement_y(
    motion_type::PeriodicRotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return sind.(γ) .* cosd.(β) .* x +
           (sind.(γ) .* sind.(β) .* sind.(α) .+ cosd.(γ) .* cosd.(α)) .* y +
           (sind.(γ) .* sind.(β) .* cosd.(α) .- cosd.(γ) .* sind.(α)) .* z .- y
end

function displacement_z(
    motion_type::PeriodicRotation{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time_triangular(t, motion_type.period, motion_type.asymmetry)
    α = t_unit .* motion_type.pitch
    β = t_unit .* motion_type.roll
    γ = t_unit .* motion_type.yaw
    return -sind.(β) .* x + cosd.(β) .* sind.(α) .* y .- z
end

function times(motion_type::PeriodicRotation)
    return [0, motion_type.period * motion_type.asymmetry, motion_type.period]
end
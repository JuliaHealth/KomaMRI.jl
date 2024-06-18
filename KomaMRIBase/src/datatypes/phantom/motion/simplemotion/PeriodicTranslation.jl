@doc raw"""
    periodic_translation = PeriodicTranslation(dx, dy, dz, period, asymmetry)

PeriodicTranslation motion struct. It produces a periodic translation of the phantom in the three directions x, y and z.
The amplitude of the oscillation will be defined by dx, dy and dz

# Arguments
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z
- `period`: (`::Real`, `[s]`) period 
- `asymmetry`: (`::Real`)  asymmetry factor, between 0 and 1

# Returns
- `periodic_translation`: (`::PeriodicTranslation`) PeriodicTranslation struct

"""
@with_kw struct PeriodicTranslation{T<:Real} <: SimpleMotionType{T}
    dx           :: T
    dy           :: T
    dz           :: T
    period::T    = typeof(dx)(0.0)
    asymmetry::T = typeof(dx)(0.5)
end

function displacement_x(
    motion_type::PeriodicTranslation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* motion_type.dx
end

function displacement_y(
    motion_type::PeriodicTranslation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* motion_type.dy
end

function displacement_z(
    motion_type::PeriodicTranslation{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* motion_type.dz
end

function times(motion_type::PeriodicTranslation)
    return [0, motion_type.period * motion_type.asymmetry, motion_type.period]
end
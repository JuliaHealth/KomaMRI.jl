@doc raw"""
    periodic_translation = PeriodicTranslation(period, asymmetry, dx, dy, dz)

PeriodicTranslation motion struct. It produces a periodic translation of the phantom in the three directions x, y and z.
The amplitude of the oscillation will be defined by dx, dy and dz


# Arguments
- `period`: (`::Real`, `[s]`) period 
- `asymmetry`: (`::Real`)  asymmetry factor, between 0 and 1
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z

# Returns
- `periodic_translation`: (`::PeriodicTranslation`) PeriodicTranslation struct


"""

@with_kw struct PeriodicTranslation{T<:Real} <: SimpleMotionType{T}
    period::T
    asymmetry::T = typeof(period)(0.5)
    dx::T        = typeof(period)(0.0)
    dy::T        = typeof(period)(0.0)
    dz::T        = typeof(period)(0.0)
end

displacement_x(motion_type::PeriodicTranslation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* motion_type.dx
end

displacement_y(motion_type::PeriodicTranslation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* motion_type.dy
end

displacement_z(motion_type::PeriodicTranslation{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* motion_type.dz
end

get_time_nodes(motion_type::PeriodicTranslation) = begin
    return [0, motion_type.period * motion_type.asymmetry, motion_type.period]
end
@doc raw"""
"""

@with_kw struct PeriodicHeartBeat{T<:Real} <: SimpleMotionType{T} 
    period::T
    asymmetry::T              = typeof(period)(0.5)
    circunferential_strain::T = typeof(period)(0.0)
    radial_strain::T          = typeof(period)(0.0)
    longitudinal_strain::T    = typeof(period)(0.0)
end 

displacement_x(motion_type::PeriodicHeartBeat{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    r = sqrt.(x.^2 + y.^2)
    θ = atan.(y, x)
    Δ_circunferential = motion_type.circunferential_strain * maximum(r)
    Δ_radial = -motion_type.radial_strain * (maximum(r) .- r)
    Δr = t_unit .* (Δ_circunferential .+ Δ_radial)
    # Map negative radius to r=0
    neg  = (r .+ Δr) .< 0
    Δr = (.!neg) .* Δr
    Δr .-= neg .* r
    return Δr .* cos.(θ)
end

displacement_y(motion_type::PeriodicHeartBeat{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    r = sqrt.(x.^2 + y.^2)
    θ = atan.(y, x)
    Δ_circunferential = motion_type.circunferential_strain * maximum(r)
    Δ_radial = -motion_type.radial_strain * (maximum(r) .- r)
    Δr = t_unit .* (Δ_circunferential .+ Δ_radial)
    # Map negative radius to r=0
    neg  = (r .+ Δr) .< 0
    Δr = (.!neg) .* Δr
    Δr .-= neg .* r
    return Δr .* sin.(θ)
end

displacement_z(motion_type::PeriodicHeartBeat{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = normalize_time_triangular(t, motion_type.period, motion_type.asymmetry)
    return t_unit .* (z .* motion_type.longitudinal_strain) 
end

get_time_nodes(motion_type::PeriodicHeartBeat) = begin
    return [0, motion_type.period * motion_type.asymmetry, motion_type.period]
end
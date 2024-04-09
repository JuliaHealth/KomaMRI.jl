@doc raw"""
"""

@with_kw struct HeartBeat{T<:Real} <: SimpleMotionType{T} 
    ti::T
    tf::T
    circunferential_strain::T = typeof(ti)(0.0)
    radial_strain::T          = typeof(ti)(0.0)
    longitudinal_strain::T    = typeof(ti)(0.0)
end 

displacement_x(motion_type::HeartBeat{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
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

displacement_y(motion_type::HeartBeat{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
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

displacement_z(motion_type::HeartBeat{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = min.(max.((t .- motion_type.ti)./(motion_type.tf - motion_type.ti), 0), 1)
    return t_unit .* (z .* motion_type.longitudinal_strain) 
end

get_time_nodes(motion_type::HeartBeat) = begin
    return [motion_type.ti, motion_type.tf]
end
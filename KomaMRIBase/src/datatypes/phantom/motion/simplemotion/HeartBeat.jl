@doc raw"""
"""

@with_kw struct HeartBeat{T<:Real} <: SimpleMotionType{T} 
    t_end::T
    t_start::T                = typeof(t_end)(0.0)
    circunferential_strain::T = typeof(t_end)(0.0)
    radial_strain::T          = typeof(t_end)(0.0)
    longitudinal_strain::T    = typeof(t_end)(0.0)
end 

is_composable(motion_type::HeartBeat) = true

displacement_x(motion_type::HeartBeat{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = motion_type.t_end == motion_type.t_start ? t .>= motion_type.t_start : min.(max.((t .- motion_type.t_start)./(motion_type.t_end - motion_type.t_start), 0), 1) 
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

displacement_y(motion_type::HeartBeat{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = motion_type.t_end == motion_type.t_start ? t .>= motion_type.t_start : min.(max.((t .- motion_type.t_start)./(motion_type.t_end - motion_type.t_start), 0), 1) 
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

displacement_z(motion_type::HeartBeat{T}, x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, t::AbstractArray{T}) where {T<:Real} = begin
    t_unit = motion_type.t_end == motion_type.t_start ? t .>= motion_type.t_start : min.(max.((t .- motion_type.t_start)./(motion_type.t_end - motion_type.t_start), 0), 1) 
    return t_unit .* (z .* motion_type.longitudinal_strain) 
end

time_nodes(motion_type::HeartBeat) = begin
    return [motion_type.t_start, motion_type.t_end]
end
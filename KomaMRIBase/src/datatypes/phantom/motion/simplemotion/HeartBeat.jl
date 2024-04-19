@doc raw"""
"""

@with_kw struct HeartBeat{T<:Real} <: SimpleMotionType{T}
    circumferential_strain :: T
    radial_strain          :: T
    longitudinal_strain::T = typeof(circumferential_strain)(0.0)
    t_start::T             = typeof(circumferential_strain)(0.0)
    t_end::T               = typeof(circumferential_strain)(0.0)
    @assert t_end >= t_start "t_end must be major or equal than t_start"
end

is_composable(motion_type::HeartBeat) = true

function displacement_x(
    motion_type::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    r = sqrt.(x .^ 2 + y .^ 2)
    θ = atan.(y, x)
    Δ_circunferential = motion_type.circumferential_strain * maximum(r)
    Δ_radial = -motion_type.radial_strain * (maximum(r) .- r)
    Δr = t_unit .* (Δ_circunferential .+ Δ_radial)
    # Map negative radius to r=0
    neg = (r .+ Δr) .< 0
    Δr = (.!neg) .* Δr
    Δr .-= neg .* r
    return Δr .* cos.(θ)
end

function displacement_y(
    motion_type::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    r = sqrt.(x .^ 2 + y .^ 2)
    θ = atan.(y, x)
    Δ_circunferential = motion_type.circumferential_strain * maximum(r)
    Δ_radial = -motion_type.radial_strain * (maximum(r) .- r)
    Δr = t_unit .* (Δ_circunferential .+ Δ_radial)
    # Map negative radius to r=0
    neg = (r .+ Δr) .< 0
    Δr = (.!neg) .* Δr
    Δr .-= neg .* r
    return Δr .* sin.(θ)
end

function displacement_z(
    motion_type::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    return t_unit .* (z .* motion_type.longitudinal_strain)
end

time_nodes(motion_type::HeartBeat) = begin
    return [motion_type.t_start, motion_type.t_end]
end
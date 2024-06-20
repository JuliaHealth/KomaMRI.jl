@doc raw"""
    heartbeat = HeartBeat(t_start, t_end, dx, dy, dz)

HeartBeat struct. It produces a heartbeat-like motion, characterised by three types of strain:
Circumferential, Radial and Longitudinal

# Arguments
- `circumferential_strain`: (`::Real`) contraction parameter
- `radial_strain`: (`::Real`) contraction parameter
- `longitudinal_strain`: (`::Real`) contraction parameter
- `t_start`: (`::Real`, `[s]`) initial time 
- `t_end`: (`::Real`, `[s]`) final time 

# Returns
- `heartbeat`: (`::HeartBeat`) HeartBeat struct

# Examples
```julia-repl
julia> hb = HeartBeat(circumferential_strain=-0.3, radial_strain=-0.2, longitudinal_strain=0.0, t_start=0.2, t_end=0.5)
```
"""
@with_kw struct HeartBeat{T<:Real} <: SimpleMotionType{T}
    circumferential_strain :: T
    radial_strain          :: T
    longitudinal_strain::T = typeof(circumferential_strain)(0.0)
    t_start::T             = typeof(circumferential_strain)(0.0)
    t_end::T               = typeof(circumferential_strain)(0.0)
    @assert t_end >= t_start "t_end must be greater or equal than t_start"
end

is_composable(motion_type::HeartBeat) = true

function displacement_x!(
    ux::AbstractArray{T},
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
    ux .= Δr .* cos.(θ)
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
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
    uy .= Δr .* sin.(θ)
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    motion_type::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion_type.t_start, motion_type.t_end)
    uz .= t_unit .* (z .* motion_type.longitudinal_strain)
    return nothing
end

times(motion_type::HeartBeat) = begin
    return [motion_type.t_start, motion_type.t_end]
end
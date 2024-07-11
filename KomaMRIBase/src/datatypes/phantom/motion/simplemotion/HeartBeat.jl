@doc raw"""
    heartbeat = HeartBeat(times, circumferential_strain, radial_strain, longitudinal_strain)

HeartBeat struct. It produces a heartbeat-like motion, characterised by three types of strain:
Circumferential, Radial and Longitudinal

# Arguments
- `times`: (`::TimeScale{T<:Real}`, `[s]`) time scale
- `circumferential_strain`: (`::Real`) contraction parameter
- `radial_strain`: (`::Real`) contraction parameter
- `longitudinal_strain`: (`::Real`) contraction parameter

# Returns
- `heartbeat`: (`::HeartBeat`) HeartBeat struct

# Examples
```julia-repl
julia> hb = HeartBeat(times=Periodic(period=1.0, asymmetry=0.3), circumferential_strain=-0.3, radial_strain=-0.2, longitudinal_strain=0.0)
```
"""
@with_kw struct HeartBeat{T<:Real, TS<:TimeScale{T}} <: SimpleMotion{T}
    times                  :: TS
    circumferential_strain :: T
    radial_strain          :: T
    longitudinal_strain    :: T = typeof(circumferential_strain)(0.0)
end

is_composable(motion::HeartBeat) = true

function displacement_x!(
    ux::AbstractArray{T},
    motion::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.times)
    r = sqrt.(x .^ 2 + y .^ 2)
    θ = atan.(y, x)
    Δ_circunferential = motion.circumferential_strain * maximum(r)
    Δ_radial = -motion.radial_strain * (maximum(r) .- r)
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
    motion::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.times)
    r = sqrt.(x .^ 2 + y .^ 2)
    θ = atan.(y, x)
    Δ_circunferential = motion.circumferential_strain * maximum(r)
    Δ_radial = -motion.radial_strain * (maximum(r) .- r)
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
    motion::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    t_unit = unit_time(t, motion.times)
    uz .= t_unit .* (z .* motion.longitudinal_strain)
    return nothing
end
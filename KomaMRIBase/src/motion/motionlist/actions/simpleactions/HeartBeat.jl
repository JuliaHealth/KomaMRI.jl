@doc raw"""
    heartbeat = HeartBeat(circumferential_strain, radial_strain, longitudinal_strain)

HeartBeat struct. It produces a heartbeat-like motion, characterised by three types of strain:
Circumferential, Radial and Longitudinal

# Arguments
- `circumferential_strain`: (`::Real`) contraction parameter
- `radial_strain`: (`::Real`) contraction parameter
- `longitudinal_strain`: (`::Real`) contraction parameter

# Returns
- `heartbeat`: (`::HeartBeat`) HeartBeat struct

# Examples
```julia-repl
julia> hb = HeartBeat(circumferential_strain=-0.3, radial_strain=-0.2, longitudinal_strain=0.0)
```
"""
@with_kw struct HeartBeat{T<:Real} <: SimpleAction{T}
    circumferential_strain :: T
    radial_strain          :: T
    longitudinal_strain    :: T
end

is_composable(action::HeartBeat) = true

function displacement_x!(
    ux::AbstractArray{T},
    action::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    r = sqrt.(x .^ 2 + y .^ 2)
    θ = atan.(y, x)
    Δ_circunferential = action.circumferential_strain * maximum(r)
    Δ_radial = -action.radial_strain * (maximum(r) .- r)
    Δr = t .* (Δ_circunferential .+ Δ_radial)
    # Map negative radius to r=0
    neg = (r .+ Δr) .< 0
    Δr = (.!neg) .* Δr
    Δr .-= neg .* r
    ux .= Δr .* cos.(θ)
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
    action::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    r = sqrt.(x .^ 2 + y .^ 2)
    θ = atan.(y, x)
    Δ_circunferential = action.circumferential_strain * maximum(r)
    Δ_radial = -action.radial_strain * (maximum(r) .- r)
    Δr = t .* (Δ_circunferential .+ Δ_radial)
    # Map negative radius to r=0
    neg = (r .+ Δr) .< 0
    Δr = (.!neg) .* Δr
    Δr .-= neg .* r
    uy .= Δr .* sin.(θ)
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    action::HeartBeat{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    uz .= t .* (z .* action.longitudinal_strain)
    return nothing
end
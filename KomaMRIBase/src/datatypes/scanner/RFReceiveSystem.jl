"""
Supertype for receive-coil sensitivity models.

Subtypes must implement [`get_n_coils`](@ref) and [`get_sens`](@ref).
"""
abstract type AbstractRFReceiveSystem end

"""Receive model with one spatially uniform channel."""
struct UniformCoilSens <: AbstractRFReceiveSystem end

"""
    ArbitraryCoilSens(x, y, z, coil_sens)

Receive model defined by complex sensitivity samples on a Cartesian grid.
`coil_sens` must have size `(length(x), length(y), length(z), ncoils)`.
Sensitivities outside the sampled grid are zero.
"""
struct ArbitraryCoilSens{
    X<:AbstractVector,
    Y<:AbstractVector,
    Z<:AbstractVector,
    S<:AbstractArray{<:Complex,4},
} <: AbstractRFReceiveSystem
    x::X
    y::Y
    z::Z
    coil_sens::S
end

"""
    BirdcageCoilSens(; ncoils=8, radius=0.20, L=0.30)

Analytic receive model with `ncoils` longitudinal birdcage rungs arranged on a
cylinder of `radius` and half-length `L`, in metres. Each rung is an ideal
finite wire, so `radius` must place the conductors outside the simulated region.
"""
Base.@kwdef struct BirdcageCoilSens <: AbstractRFReceiveSystem
    ncoils::Int = 8
    radius::Float64 = 0.20
    L::Float64 = 0.30
end

"""Return the number of receive channels in `receiver`."""
get_n_coils(::UniformCoilSens) = 1
get_n_coils(receiver::BirdcageCoilSens) = receiver.ncoils
get_n_coils(receiver::ArbitraryCoilSens) = size(receiver.coil_sens, 4)

@inline function birdcage_sensitivity(x, y, z, coil, radius, L, ncoils)
    ϕ = oftype(x, 2π) * coil / ncoils
    xn = radius * cos(ϕ)
    yn = radius * sin(ϕ)
    rn = sqrt((x - xn)^2 + (y - yn)^2 + eps(x))
    Bϕ = (1 / rn) * (
        (L - z) / sqrt(rn^2 + (L - z)^2) +
        (L + z) / sqrt(rn^2 + (L + z)^2)
    )
    return complex(-Bϕ * (y - yn) / rn, -Bϕ * (x - xn) / rn)
end

"""
    get_sens(receiver, x, y, z)

Evaluate the complex receive sensitivities at the supplied spin positions.
The result has size `(length(x), get_n_coils(receiver))`.
"""
function get_sens(receiver::BirdcageCoilSens, x, y, z)
    T = eltype(x)
    radius = T(receiver.radius)
    L = T(receiver.L)
    ncoils = get_n_coils(receiver)
    nspins = length(x)
    sens = similar(x, Complex{T}, nspins, ncoils)
    coils = similar(x, T, ncoils)
    coils .= T.(0:(ncoils - 1))
    sens .= birdcage_sensitivity.(
        reshape(x, :, 1),
        reshape(y, :, 1),
        reshape(z, :, 1),
        reshape(coils, 1, :),
        radius,
        L,
        T(ncoils),
    )
    return sens
end

function get_sens(receiver::ArbitraryCoilSens, x, y, z)
    ncoils = get_n_coils(receiver)
    sens = similar(x, eltype(receiver.coil_sens), length(x), ncoils)
    for coil in axes(receiver.coil_sens, 4)
        interpolation = linear_interpolation(
            (receiver.x, receiver.y, receiver.z),
            @view(receiver.coil_sens[:, :, :, coil]);
            extrapolation_bc=zero(eltype(receiver.coil_sens)),
        )
        sens[:, coil] .= interpolation.(x, y, z)
    end
    return sens
end

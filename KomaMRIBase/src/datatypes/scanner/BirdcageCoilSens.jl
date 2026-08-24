"""
    BirdcageCoilSens(; ncoils=8, radius=0.20, L=0.30)

Analytic receive model with `ncoils` longitudinal birdcage rungs arranged on a
cylinder of `radius` and half-length `L`, in metres. Each rung is an ideal
finite wire, so `radius` must place the conductors outside the simulated region.
"""
Base.@kwdef struct BirdcageCoilSens{R<:Real,LT<:Real} <: AbstractRFReceiveSystem
    ncoils::Int = 8
    radius::R = 0.20
    L::LT = 0.30
end

get_n_coils(receiver::BirdcageCoilSens) = receiver.ncoils

# Scalar field equations shared by CPU broadcasting and GPU kernel evaluation.
@inline function birdcage_sensitivity(x, y, z, coil, radius, L, ncoils)
    ϕ = oftype(x, 2π) * coil / ncoils
    xn = radius * cos(ϕ)
    yn = radius * sin(ϕ)
    return birdcage_wire_sensitivity(x, y, z, xn, yn, L)
end

@inline function birdcage_wire_sensitivity(x, y, z, xn, yn, L)
    rn = sqrt((x - xn)^2 + (y - yn)^2 + eps(x))
    Bϕ = (1 / rn) * (
        (L - z) / sqrt(rn^2 + (L - z)^2) +
        (L + z) / sqrt(rn^2 + (L + z)^2)
    )
    return complex(-Bϕ * (y - yn) / rn, -Bϕ * (x - xn) / rn)
end

function get_sens(receiver::BirdcageCoilSens, x, y, z)
    T = eltype(x)
    sens = similar(x, Complex{T}, length(x), get_n_coils(receiver))
    return get_sens!(sens, receiver, x, y, z)
end

function get_sens!(sens, receiver::BirdcageCoilSens, x, y, z)
    T = eltype(x)
    radius = T(receiver.radius)
    L = T(receiver.L)
    ncoils = get_n_coils(receiver)
    sens .= birdcage_sensitivity.(
        reshape(x, :, 1),
        reshape(y, :, 1),
        reshape(z, :, 1),
        reshape(0:(ncoils - 1), 1, :),
        radius,
        L,
        T(ncoils),
    )
    return sens
end

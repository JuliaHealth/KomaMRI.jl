abstract type AbstractRFReceiveSystem end

struct UniformCoilSens <: AbstractRFReceiveSystem end

struct ArbitraryCoilSens{
    C<:AbstractArray,
    S<:AbstractArray{<:Complex},
} <: AbstractRFReceiveSystem
    x::C
    y::C
    z::C
    coil_sens::S
end

Base.@kwdef struct BirdcageCoilSens <: AbstractRFReceiveSystem
    ncoils::Int = 8
    radius::Float64 = 0.20
    L::Float64 = 0.30
end

get_n_coils(::AbstractRFReceiveSystem) = 1
get_n_coils(receiver::BirdcageCoilSens) = receiver.ncoils
get_n_coils(receiver::ArbitraryCoilSens) = size(receiver.coil_sens, 4)

function get_sens(receiver::BirdcageCoilSens, x, y, z)
    T = eltype(x)
    radius = T(receiver.radius)
    L = T(receiver.L)
    ncoils = get_n_coils(receiver)
    nspins = length(x)
    sens = similar(x, Complex{T}, nspins, ncoils)
    ϵ = eps(T)

    for n in 1:ncoils
        ϕn = T(2π) * (n - 1) / ncoils
        xn = radius * cos(ϕn)
        yn = radius * sin(ϕn)
        rn = @. sqrt((x - xn)^2 + (y - yn)^2 + ϵ)
        Bϕ = @. (1 / rn) * ((L - z) / sqrt(rn^2 + (L - z)^2) + (L + z) / sqrt(rn^2 + (L + z)^2))
        Bx = @. -Bϕ * (y - yn) / rn
        By = @. Bϕ * (x - xn) / rn
        sens[:, n] .= Bx - im * By
    end
    return sens
end

function get_sens(receiver::ArbitraryCoilSens, x, y, z)
    sens = similar(receiver.coil_sens, length(x), get_n_coils(receiver))
    for coil in axes(receiver.coil_sens, 4)
        base_itp = GriddedInterpolation(
            (receiver.x, receiver.y, receiver.z),
            receiver.coil_sens[:, :, :, coil],
            Gridded(Linear()),
        )
        itp = extrapolate(base_itp, zero(eltype(receiver.coil_sens)))
        sens[:, coil] .= itp.(x, y, z)
    end
    return sens
end

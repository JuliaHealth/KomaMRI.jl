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

function sensitivity_interpolator(receiver::ArbitraryCoilSens)
    T = eltype(receiver.x)
    coils = similar(receiver.x, T, get_n_coils(receiver))
    coils .= axes(receiver.coil_sens, 4)
    base_itp = GriddedInterpolation(
        (receiver.x, receiver.y, receiver.z, coils),
        receiver.coil_sens,
        Gridded(Linear()),
    )
    return extrapolate(base_itp, zero(eltype(receiver.coil_sens)))
end

# Evaluate every coil at the supplied spin positions.
function interpolate_sensitivities(itp, x, y, z)
    sens = similar(x, eltype(itp), length(x), size(itp, 4))
    sens .= itp.(
        reshape(x, :, 1),
        reshape(y, :, 1),
        reshape(z, :, 1),
        reshape(last(itp.itp.knots), 1, :),
    )
    return sens
end

get_sens(receiver::ArbitraryCoilSens, x, y, z) =
    interpolate_sensitivities(sensitivity_interpolator(receiver), x, y, z)

abstract type AbstractReceiveSystem end

struct UniformCoilSens <: AbstractReceiveSystem end

const ComplexScalarOrArray{T} = Union{Complex{T}, AbstractArray{Complex{T}}} where {T}
struct ArbitraryRFRxCoils{T,
                        T1 <: ComplexScalarOrArray{T},
                        T2 <: ComplexScalarOrArray{T}
                        } <: AbstractReceiveSystem 
    x::AbstractVector{T}
    y::AbstractVector{T}
    z::AbstractVector{T}
    coil_sens::T1  
    B1⁺::T2
end

Base.@kwdef struct BirdcageCoilSens <: AbstractReceiveSystem
    ncoils::Int = 8
    radius::Float64 = 0.20
    L::Float64 = 0.30
end

get_n_coils(::AbstractReceiveSystem) = 1
get_n_coils(receiver::BirdcageCoilSens) = receiver.ncoils
get_n_coils(receiver::ArbitraryRFRxCoils) = size(receiver.coil_sens, 4)

function get_sens(receiver::BirdcageCoilSens, x, y, z)
    T = eltype(x)
    radius = T(receiver.radius)
    L = T(receiver.L)
    ncoils = get_n_coils(receiver)
    nspins = length(x)
    sens = Matrix{Complex{T}}(undef, nspins, ncoils)
    ϵ = eps(T)

    for n in 1:ncoils
        ϕn = T(2π) * (n - 1) / ncoils
        xn = radius * cos(ϕn)
        yn = radius * sin(ϕn)
        rn = @. sqrt((x - xn)^2 + (y - yn)^2 + ϵ)
        Bϕ = @. (1 / rn) * ((L - z) / sqrt(rn^2 + (L - z)^2) + (L + z) / sqrt(rn^2 + (L + z)^2))
        Bx = @. -Bϕ * (y - yn) / rn
        By = @. Bϕ * (x - xn) / rn
        sens[:, n] .=  Bx - im * By
    end
    return sens
end

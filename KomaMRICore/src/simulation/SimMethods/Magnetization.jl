"""
    mag = Mag(xy::Complex, z::Real)

The Magnetization struct.

# Arguments
- `xy`: (`::Complex{Float64}`) magnetization of a spin in the xy plane
- `z`: (`::Real`) magnetization of a spin in the z plane

# Returns
- `mag`: (`::Mag`) Magnetization struct
"""
mutable struct Mag{XYType<:AbstractVector,ZType<:AbstractVector} <: SpinStateRepresentation
    xy::XYType
    z::ZType
end

# Required indexing operations
# M[i]
Base.getindex(M::Mag, i::Integer) = Mag(M.xy[i,:], M.z[i,:])
# M[a:b]
Base.getindex(M::Mag, i) = Mag(M.xy[i], M.z[i])
Base.view(M::Mag, i) = @views Mag(M.xy[i], M.z[i])

@inline function spinor_half_angle(θ2::T) where {T<:Real}
    if θ2 <= sqrt(eps(T))
        θ4 = θ2 * θ2
        return T(1) - θ2 / T(8) + θ4 / T(384),
               T(0.5) - θ2 / T(48) + θ4 / T(3840)
    end
    θ = sqrt(θ2)
    s, c = sincos(θ / T(2))
    return c, s / θ
end

function set_rotation_spinor!(α, β, θxy, θz)
    @inbounds for i in eachindex(α, β, θxy, θz)
        c, scale = spinor_half_angle(abs2(θxy[i]) + θz[i]^2)
        α[i] = complex(c, -θz[i] * scale)
        β[i] = complex(imag(θxy[i]) * scale, -real(θxy[i]) * scale)
    end
    return nothing
end

calc_mag_norm!(norm, M::Mag) = nothing
function calc_mag_norm!(norm::AbstractVector{Float32}, M::Mag)
    @. norm = sqrt(abs2(M.xy) + M.z^2)
    return nothing
end

restore_mag_norm!(norm, M::Mag) = nothing
function restore_mag_norm!(norm::AbstractVector{Float32}, M::Mag)
    @. norm = norm / sqrt(abs2(M.xy) + M.z^2)
    @. M.xy = M.xy * norm
    @. M.z = M.z * norm
    return nothing
end

# Definition of rotation Spinor×SpinStateRepresentation
@doc raw"""
Spinor (\alpha, \beta) × Magnetization (Mx + i My, Mz)

A vector M = (Mx,My,Mz) can be expressed as a complex 2x2 matrix

M = [Mz Mxy⋆;

------- Mxy  -Mz],	with Mxy = Mx + i My.

Then, to operate with a Spinor M+=RMR⋆, or (α,β)×(Mxy,Mz) = (Mxy+,Mz+) with

Mxy+ = 2α⋆βMz+(α⋆)²Mxy-β²Mxy⋆

and

Mz+ = (|α|² - |β|²)Mz-α⋆ β⋆ Mxy-αβMxy⋆ .

Pauly, J., Le Roux, P., Nishimura, D., & Macovski, A. (1991).
Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm
(NMR imaging).
IEEE Transactions on Medical Imaging, 10(1), 53-65. doi:10.1109/42.75611
"""
mul!(s::Spinor, M::Mag) = begin
    M_aux = Mag(
        2 .*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs2.(s.α) .- abs2.(s.β)).*M.z .- 2 .*real.(s.α.*s.β.*conj.(M.xy))
     )
    M.xy .= M_aux.xy
    M.z  .= M_aux.z
end
@inline mul!(s::Spinor, M::Mag, Maux_xy, Maux_z) = begin
    @. Maux_xy = 2*conj(s.α)*s.β*M.z+conj(s.α)^2*M.xy-s.β^2*conj(M.xy)
    @. Maux_z = (abs2(s.α) - abs2(s.β))*M.z - 2*real(s.α*s.β*conj(M.xy))
    @. M.xy = Maux_xy
    @. M.z = Maux_z
end
@inline *(s::Spinor, M::Mag) = begin
    Mag(
        2 .*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs2.(s.α) .- abs2.(s.β)).*M.z .- 2 .*real.(s.α.*s.β.*conj.(M.xy))
     )
end

"""
    mag = Mag(xy::Complex, z::Real)

The Magnetization struct.

# Arguments
- `xy`: (`::Complex{Float64}`) magnetization of a spin in the xy plane
- `z`: (`::Real`) magnetization of a spin in the z plane

# Returns
- `mag`: (`::Mag`) Magnetization struct
"""
mutable struct Mag{T<:Real} <: SpinStateRepresentation{T}
    xy::AbstractVector{Complex{T}}
    z::AbstractVector{T}
end

# Required indexing operations
# M[i]
Base.getindex(M::Mag, i::Integer) = Mag(M.xy[i,:], M.z[i,:])
# M[a:b]
Base.getindex(M::Mag, i) = Mag(M.xy[i], M.z[i])
Base.view(M::Mag, i) = @views Mag(M.xy[i], M.z[i])

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
mul!(s::Spinor{T}, M::Mag) where {T<:Real} = begin
    M_aux = Mag(
        T(2) .*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-T(2) .*real.(s.α.*s.β.*conj.(M.xy))
     )
    M.xy .= M_aux.xy
    M.z  .= M_aux.z
end
mul!(s::Spinor{T}, M::Mag, Maux_xy, Maux_z) where {T<:Real} = begin
    @. Maux_xy = T(2)*conj(s.α)*s.β*M.z+conj(s.α)^2*M.xy-s.β^2*conj(M.xy)
    @. Maux_z = (abs(s.α)^2 -abs(s.β)^2)*M.z-T(2) *real(s.α*s.β*conj(M.xy))
    @. M.xy = Maux_xy
    @. M.z = Maux_z
end
*(s::Spinor{T}, M::Mag) where {T<:Real} = begin
    Mag(
        T(2) .*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-T(2) .*real.(s.α.*s.β.*conj.(M.xy))
     )
end

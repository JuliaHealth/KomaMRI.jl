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
    xy::AbstractVector
    z::AbstractVector
end

# Convenience constructor used by Functors/concrete callers; Reactant paths use Mag{T}(...) directly
Mag(xy::AbstractVector{Complex{T}}, z::AbstractVector{T}) where {T<:Real} = Mag{T}(xy, z)

# Required indexing operations
# M[i]
Base.getindex(M::Mag{T}, i::Integer) where {T} = Mag{T}(M.xy[i,:], M.z[i,:])
# M[a:b]
Base.getindex(M::Mag{T}, i) where {T} = Mag{T}(M.xy[i], M.z[i])
Base.view(M::Mag{T}, i) where {T} = @views Mag{T}(M.xy[i], M.z[i])

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
mul!(s::Spinor, M::Mag{T}) where {T} = begin
    M_aux = Mag{T}(
        2 .*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-2 .*real.(s.α.*s.β.*conj.(M.xy))
     )
    M.xy .= M_aux.xy
    M.z  .= M_aux.z
end
mul!(s::Spinor, M::Mag, Maux_xy, Maux_z) = begin
    @. Maux_xy = 2*conj(s.α)*s.β*M.z+conj(s.α)^2*M.xy-s.β^2*conj(M.xy)
    @. Maux_z = (abs(s.α)^2 -abs(s.β)^2)*M.z-2 *real(s.α*s.β*conj(M.xy))
    @. M.xy = Maux_xy
    @. M.z = Maux_z
end
*(s::Spinor, M::Mag{T}) where {T} = begin
    Mag{T}(
        2 .*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-2 .*real.(s.α.*s.β.*conj.(M.xy))
     )
end

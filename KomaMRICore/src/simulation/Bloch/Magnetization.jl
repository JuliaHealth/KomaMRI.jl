"""
    mag = Mag(xy, z)

The Magnetization struct.

# Arguments
- `xy`: (`::AbstractVector{Complex{Real}}`) magnetization of spins in the xy plane
- `z`: (`::AbstractVector{Real}`) magnetization of spins in the z plane

# Returns
- `mag`: (`::SpinStateRepresentation{Real}`) Magnetization struct
"""
mutable struct Mag{T<:Real} <: SpinStateRepresentation{T}
    xy::AbstractVector{Complex{T}}
    z::AbstractVector{T}
end

# Required indexing operations
# M[i]
Base.getindex(M::Mag, i::Integer) = Mag(M.xy[i,:], M.z[i,:])
# M[a:b]
Base.getindex(M::Mag, i::UnitRange) = Mag(M.xy[i], M.z[i])
Base.view(M::Mag, i::UnitRange) = @views Mag(M.xy[i], M.z[i])

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
        2*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-2*real.(s.α.*s.β.*conj.(M.xy))
     )
    M.xy .= M_aux.xy
    M.z  .= M_aux.z
end
*(s::Spinor, M::Mag) = begin
    Mag(
        2*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-2*real.(s.α.*s.β.*conj.(M.xy))
     )
end

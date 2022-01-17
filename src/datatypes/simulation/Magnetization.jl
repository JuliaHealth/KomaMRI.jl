mutable struct Mag
    xy::Complex
    z::Real
end 
Base.show(io::IO,M::Mag) = print(io, "Mag(xy = ",round(M.xy,digits=2),", z = ",round(M.z,digits=2),")")
# Contructor
Mag(p::Phantom, dir::Symbol) = dir==:x ? Mag.(p.ρ,0) : Mag.(0,p.ρ)
# Arithmetic operations
+(M1::Mag, M2::Mag) = Mag(M1.xy + M2.xy, M1.z + M2.z) #Vector sum
*(α::Float64, M::Mag) = Mag(α*M.xy, α*M.z)
# Other operations
angle(M::Mag) = angle(M.xy)
abs(M::Mag) = abs(M.xy)
getproperty(x::Vector{Mag}, f::Symbol) = getproperty.(x,f)
# Rotation
@doc raw"""
Spinor (\alpha, \beta) × Magnetization (Mx + i My, Mz)

A vector M = (Mx,My,Mz) can be expressed as a complex 2x2 matrix

M = [Mz Mxy⋆;

------- Mxy  -Mz],	with Mxy = Mx + i My.

Then, to operate with a Spinor M+=RMR⋆, or (α,β)×(Mxy,Mz) = (Mxy+,Mz+) with

Mxy+ = 2α⋆βMz+(α⋆)²Mxy-β²Mxy⋆

and

Mz+ = (|α|² - |β|²)Mz-α⋆ β⋆ Mxy-αβMxy⋆ .
"""
*(s::Spinor, M::Mag) = begin
    # M_mat = [M.z conj(M.xy); M.xy -M.z]
    # Q = [s.α -conj(s.β); s.β conj(s.α)]
    # M⁺ = Q * M_mat * Q'
    # Mag(M⁺[2,1], real(M⁺[1,1]))
	Mag(
        2*conj(s.α)*s.β*M.z+conj(s.α)^2*M.xy-s.β^2*conj(M.xy),
        (abs(s.α)^2-abs(s.β)^2)*M.z-conj(s.α)*conj(s.β)*M.xy-s.α*s.β*conj(M.xy)
        )
end
#Operation on vector
# M0 = Mag(0,1)
# Mf = Ry(π/2)*M0
# Mf.xy, Mf.z

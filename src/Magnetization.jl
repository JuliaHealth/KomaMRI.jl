mutable struct Mag
    xy::Complex
    z::Real
end 
Base.show(io::IO,M::Mag) = print(io, "Mag(xy = ",round(M.xy,digits=2),", z = ",round(M.z,digits=2),")")
# Contructor
Mag(p::Phantom) = Mag.(0,p.ρ)
# Arithmetic operations
+(M1::Mag, M2::Mag) = Mag(M1.xy + M2.xy, M1.z + M2.z) #Vector sum
*(α::Float64, M::Mag) = Mag(α*M.xy, M.z)
angle(M::Mag) = angle(M.xy)
abs(M::Mag) = abs(M.xy)
getproperty(x::Vector{Mag}, f::Symbol) = getproperty.(x,f)
"""
Spinor × Magnetization (Mx + i My, Mz)

A vector V = (x,y,z) can be expressed as a complex 2x2 matrix

V = [z X⋆;

------- X  -z],	with X = x + i y.

Then, to operate with a Spinor V+=RVR⋆, or (α,β)×(X,z) = (X+,z+) with

X+ = 2α⋆βz+(α⋆)²X-β²X⋆

and

z+ = (|α|² - |β|²)z-α⋆ β⋆ X-αβX⋆ .
"""
*(s::Spinor, M::Mag) = begin
	Mag(
        2*conj(s.α)*s.β*M.z+conj(s.α)^2*M.xy-s.β^2*conj(M.xy),
        (abs(s.α)^2-abs(s.β)^2)*M.z-conj(s.α)*conj(s.β)*M.xy-s.α*s.β*conj(M.xy)
        )
end

# M0 = Mag(0,1)
# Mf = Ry(π/2)*M0
# Mf.xy, Mf.z

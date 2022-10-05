"""
    mag = Mag(xy::Complex, z::Real)

The Magnetization struct.

# Arguments
- `xy`: (`::Complex{Int64}`) the magnetization of a spin in the xy plane
- `z`: (`::Real`) the magnetization of a spin in the z plane

# Returns
- `mag`: (`::Mag`) the Magnetization struct
"""
mutable struct Mag
    xy::Complex
    z::Real
end

"""
    mag = Mag(p::Phantom, dir::Symbol)

Generates a Vector of Mag structs with the information of the `p` panthom struct (the proton
density) in the axis given by the `dir` symbol.

# Arguments
- `p`: (`::Phantom`) the phantom struct
- `dir`: (`::Symbol`, opts: [`:x`, `:z`]) the symbol that represents the axis of the
    magnetization

# Returns
- `mag`: (`Vector{Mag}`) the vector of Magnetization structs

# Examples
```julia-repl
julia> obj = Phantom(x=zeros(5));

julia> mag = Mag(obj, :z)
5-element Vector{Mag}:
 Mag(xy = 0.0 + 0.0im, z = 1.0)
 Mag(xy = 0.0 + 0.0im, z = 1.0)
 Mag(xy = 0.0 + 0.0im, z = 1.0)
 Mag(xy = 0.0 + 0.0im, z = 1.0)
 Mag(xy = 0.0 + 0.0im, z = 1.0)
```
"""
Mag(p::Phantom, dir::Symbol) = dir==:x ? Mag.(p.ρ, 0) : Mag.(0, p.ρ)

"""
    str = show(io::IO, x::Mag)

Displays information about the Mag struct `x` in the julia REPL.

# Arguments
- `x`: (`::Mag`) the Magnetization struct

# Returns
- `str` (`::String`) the output string message
"""
Base.show(io::IO, M::Mag) = begin
    print(io, "Mag(xy = ", round(M.xy, digits=2), ", z = ", round(M.z, digits=2), ")")
end

"""
    y = getproperty(x::Vector{Mag}, f::Symbol)

Overchages Base.getproperty(). It is meant to access properties of the Mag vector `x`
directly without needing to iterate elementwise.

# Arguments
- `x`: (`::Vector{Mag}`) the vector of Mag structs
- `f`: (`::Symbol`, opts: [`:xy`, `:z`]) the symbol that represents a property of a Mag
    struct

# Returns
- `y`: (`::Vector{Any}`) the vector with the property defined by the `f` symbol for all
    elements of the Mag vector `x`
"""
getproperty(x::Vector{Mag}, f::Symbol) = getproperty.(x, f)


# Arithmetic operations
+(M1::Mag, M2::Mag) = Mag(M1.xy + M2.xy, M1.z + M2.z) #Vector sum
*(α::Float64, M::Mag) = Mag(α*M.xy, α*M.z) #This was defined to easily define the average M_mean = sum 1/N* M
*(M::Mag, α::Float64) = Mag(α*M.xy, α*M.z)

# Other operations
angle(M::Mag) = angle(M.xy)
abs(M::Mag) = abs(M.xy)

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
    # Matricial representation:
    # M_mat = [M.z conj(M.xy); M.xy -M.z]
    # Q = [s.α -conj(s.β); s.β conj(s.α)]
    # M⁺ = Q * M_mat * Q'
    # Mag(M⁺[2,1], real(M⁺[1,1]))
	Mag(
        2*conj(s.α)*s.β*M.z+conj(s.α)^2*M.xy-s.β^2*conj(M.xy),
        (abs(s.α)^2-abs(s.β)^2)*M.z-2real(s.α*s.β*conj(M.xy))#2real(s.α*s.β*conj(M.xy))=conj(s.α)*conj(s.β)*M.xy+s.α*s.β*conj(M.xy)
        )
end

#Operation on vector
# M0 = Mag(0,1)
# Mf = Ry(π/2)*M0
# Mf.xy, Mf.z

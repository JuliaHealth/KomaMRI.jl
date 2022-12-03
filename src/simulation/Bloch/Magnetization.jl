"""
    mag = Mag(xy::Complex, z::Real)

The Magnetization struct.

# Arguments
- `xy`: (`::Complex{Float64}`) magnetization of a spin in the xy plane
- `z`: (`::Real`) magnetization of a spin in the z plane

# Returns
- `mag`: (`::Mag`) Magnetization struct
"""
mutable struct Mag{T<:Real} <: SpinStateRepresentation
    xy::AbstractVector{Complex{T}}
    z::AbstractVector{T}
end
Mag(xy::Complex{T}, z::T) where {T<:Real} = Mag([xy], [z])
Mag(xy::T, z::T) where {T<:Real} = Mag([complex(xy)], [z])

"""
    str = show(io::IO, x::Mag)

Displays information about the Mag struct `x` in the julia REPL.

# Arguments
- `x`: (`::Mag`) Magnetization struct

# Returns
- `str` (`::String`) output string message
"""
Base.show(io::IO, M::Mag) = begin
    print(io, "Mag(xy = ", round.(M.xy, digits=2), ", z = ", round.(M.z, digits=2), ")")
end

"""
    y = getproperty(x::Vector{Mag}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the Mag vector `x`
directly without needing to iterate elementwise.

# Arguments
- `x`: (`::Vector{Mag}`) vector of Mag structs
- `f`: (`::Symbol`, opts: [`:xy`, `:z`]) symbol that represents a property of a Mag
    struct

# Returns
- `y`: (`::Vector{Any}`) vector with the property defined by the `f` symbol for all
    elements of the Mag vector `x`
"""
getproperty(x::Vector{Mag}, f::Symbol) = getproperty.(x, f)

# Arithmetic operations
+(M1::Mag, M2::Mag) = Mag(M1.xy .+ M2.xy, M1.z .+ M2.z) #Vector sum
*(α::Float64, M::Mag) = Mag(α*M.xy, α*M.z) #This was defined to easily define the average M_mean = sum 1/N* M
*(M::Mag, α::Float64) = Mag(α*M.xy, α*M.z)

# Other operations
angle(M::Mag) = angle.(M.xy)
abs(M::Mag) = abs.(M.xy)
zero(::Type{Mag}) = Mag([0],[0])
Base.getindex(M::Mag, i) = Mag(M.xy[i], M.z[i])
Base.lastindex(M::Mag) = length(M.xy)
Base.view(M::Mag, i) = Mag(view(M.xy, i), view(M.z, i))
Base.length(M::Mag) = length(M.xy)
Base.iterate(M::Mag) = (Mag(M.xy[1], M.z[1]), 2)
Base.iterate(M::Mag, i::Integer) = (i <= length(M)) ? (Mag(M.xy[i], M.z[i]), i+1) : nothing

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
        2*conj.(s.α).*s.β.*M.z.+conj.(s.α).^2 .* M.xy.-s.β.^2 .*conj.(M.xy),
        (abs.(s.α).^2 .-abs.(s.β).^2).*M.z.-2*real.(s.α.*s.β.*conj.(M.xy)) #2real(s.α*s.β*conj(M.xy))=conj(s.α)*conj(s.β)*M.xy+s.α*s.β*conj(M.xy)
        )
end

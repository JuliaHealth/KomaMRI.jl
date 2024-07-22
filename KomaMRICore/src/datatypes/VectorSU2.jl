mutable struct VectorSU2{T}
    xy::AbstractVector{Complex{T}}
    z::AbstractVector{T}
end

# Required indexing operations
Base.view(V::VectorSU2, i::UnitRange) = @views VectorSU2(V.xy[i], V.z[i])

# Maths
Base.abs(V::VectorSU2) = sqrt.(abs.(V.xy).^2 + (V.z).^2)
Base.:/(V::VectorSU2, α) = VectorSU2(V.xy/α, V.z/α)
Base.:+(V1::VectorSU2, V2::VectorSU2) = VectorSU2(V1.xy .+ V2.xy, V1.z .+ V2.z)

# Definition of Spinor from VectorSU2
function calculateRot!(pre::PreallocResult{T}, Δt::T) where {T}
    # Beff = pre.B_new
    Beff = (pre.B_old + pre.B_new) / 2
    absBeff = abs(Beff)
    # absBeff[absBeff .== zero(T)] .= eps(T) # not sure why this is needed, in RF block this should be non-zero
    @. pre.φ = T(-2π * γ) * absBeff * Δt 
    @. pre.Rot.α = cos(pre.φ/2) - Complex{T}(im) * (Beff.z / absBeff) * sin(pre.φ/2)
    @. pre.Rot.β = -Complex{T}(im) * (Beff.xy / absBeff) * sin(pre.φ/2)
    return nothing
end
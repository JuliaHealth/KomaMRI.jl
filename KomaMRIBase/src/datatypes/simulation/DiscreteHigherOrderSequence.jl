"""
    seqd = DiscreteHigherOrderSequence{T, order}(G, B1, Δf, ADC, t, Δt)

A sampled version of a higher order Sequence struct, containing vectors for event amplitudes at specified
times. DiscreteHigherOrderSequence is the struct used for simulation.

# Arguments
- `G`: (`::AbstractMatrix{T<:Real}`) Spatial harmonic gradient matrix, with each row corresponding to a gradient order
    Order -1 Gx, Gy, Gz
    Order 0  B0
    Order 1  B0 Gx, Gy, Gz
    Order 2  B0 Gx, Gy, Gz, Gxy, Gyz, G(2z^2 - x^2 - y^2), Gxz, G(x^2 - y^2)
- `B1`: (`::AbstractVector{Complex{T<:Real}}`, `[T]`) RF amplitude vector
- `Δf`: (`::AbstractVector{T<:Real}`, `[Hz]`) RF carrier frequency displacement vector
- `ADC`: (`::AbstractVector{Bool}`) ADC sample vector
- `t`: (`::AbstractVector{T<:Real}`, `[s]`) time vector
- `Δt`: (`::AbstractVector{T<:Real}`, `[s]`) delta time vector

# Returns
- `seqd`: (`::DiscreteHigherOrderSequence`) DiscreteHigherOrderSequence struct
"""
struct DiscreteHigherOrderSequence{T<:Real, order<:Integer} <: AbstractDiscreteSequence{T}
    G::AbstractMatrix{T}
    B1::AbstractVector{Complex{T}}
    Δf::AbstractVector{T}
    ADC::AbstractVector{Bool}
    t::AbstractVector{T}
    Δt::AbstractVector{T}
end

Base.getindex(seq::DiscreteHigherOrderSequence{T, order}, i::Integer) where {T<:Real, order<:Integer} = begin
    DiscreteHigherOrderSequence{T, order}(seq.G[:, i],
                                          seq.B1[i, :],
                                          seq.Δf[i, :],
                                          seq.ADC[i, :],
                                          seq.t[i, :],
                                          seq.Δt[i, :])
end
Base.getindex(seq::DiscreteHigherOrderSequence{T, order}, i::UnitRange) where {T<:Real, order<:Integer} = begin
    DiscreteHigherOrderSequence{T, order}(seq.G[:, i],
                                          seq.B1[i],
                                          seq.Δf[i],
                                          seq.ADC[i],
                                          seq.t[i],
                                          seq.Δt[i.start:i.stop-1])
end
Base.view(seq::DiscreteHigherOrderSequence{T, order}, i::UnitRange) where {T<:Real, order<:Integer} = begin
    @views DiscreteHigherOrderSequence{T, order}(seq.G[:, i],
                                                 seq.B1[i],
                                                 seq.Δf[i],
                                                 seq.ADC[i],
                                                 seq.t[i],
                                                 seq.Δt[i.start:i.stop-1])
end

is_GR_on(seq::DiscreteHigherOrderSequence) =  sum(abs.(seq.G[:, 1:end-1])) != 0


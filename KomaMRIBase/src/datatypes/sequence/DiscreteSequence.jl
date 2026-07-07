# ===========================================================================
# 0. Sampled Sequence Storage
# ===========================================================================
#
# DiscreteSequence stores RF, gradient, frequency, and ADC event waveforms after a
# sequence has been sampled on a time grid.

# -- 0.1. Type definition ----------------------------------------------------
"""
    seqd = DiscreteSequence(Gx, Gy, Gz, B1, Δf, ψ, ADC, excitation_bool, t, Δt)

A sampled version of a Sequence struct, containing event amplitudes at specified times.

# Arguments
- `Gx`: (`::AbstractVector{T<:Real}`, `[T/m]`) x-gradient vector
- `Gy`: (`::AbstractVector{T<:Real}`, `[T/m]`) y-gradient vector
- `Gz`: (`::AbstractVector{T<:Real}`, `[T/m]`) z-gradient vector
- `B1`: (`::AbstractVector{Complex{T<:Real}}`, `[T]`) RF amplitude vector
- `Δf`: (`::AbstractVector{T<:Real}`, `[Hz]`) RF carrier frequency displacement vector
- `ψ`: (`::AbstractVector{T<:Real}`, `[rad]`) RF rotating-frame phase vector
- `ADC`: (`::AbstractVector{Bool}`) ADC sample vector
- `excitation_bool`: (`::AbstractVector{Bool}`) RF excitation interval vector
- `t`: (`::AbstractVector{T<:Real}`, `[s]`) time vector
- `Δt`: (`::AbstractVector{T<:Real}`, `[s]`) delta time vector

# Returns
- `seqd`: (`::DiscreteSequence`) DiscreteSequence struct
"""
struct DiscreteSequence{
    GxType<:AbstractVector,
    GyType<:AbstractVector,
    GzType<:AbstractVector,
    B1Type<:AbstractVector,
    ΔfType<:AbstractVector,
    ψType<:AbstractVector,
    ADCType<:AbstractVector{Bool},
    ExcitationType<:AbstractVector{Bool},
    tType<:AbstractVector,
}
    Gx::GxType
    Gy::GyType
    Gz::GzType
    B1::B1Type
    Δf::ΔfType
    ψ::ψType
    ADC::ADCType
    excitation_bool::ExcitationType
    t::tType
    Δt::tType
end

# -- 0.2. Storage helpers ----------------------------------------------------
function DiscreteSequence(t::AbstractVector{T}=Float64[]) where {T<:Real}
    n = length(t)
    return DiscreteSequence(
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(Complex{T}, n),
        zeros(T, n),
        zeros(T, n),
        fill(false, n),
        Bool[],
        t,
        similar(t, 0),
    )
end

table_columns(seqd) = (seqd.Gx, seqd.Gy, seqd.Gz, seqd.B1, seqd.Δf, seqd.ψ, seqd.ADC)

# -- 0.3. Indexing and iteration --------------------------------------------
Base.length(seq::DiscreteSequence) = length(seq.Δt)
Base.getindex(seq::DiscreteSequence, i::Integer) = begin
    DiscreteSequence(seq.Gx[i, :],
                     seq.Gy[i, :],
                     seq.Gz[i, :],
                     seq.B1[i, :],
                     seq.Δf[i, :],
                     seq.ψ[i, :],
                     seq.ADC[i, :],
                     seq.excitation_bool[i, :],
                     seq.t[i, :],
                     seq.Δt[i, :])
end

Base.getindex(seq::DiscreteSequence, i::UnitRange) = begin
    intervals = i.start:i.stop-1
    DiscreteSequence(seq.Gx[i],
                     seq.Gy[i],
                     seq.Gz[i],
                     seq.B1[i],
                     seq.Δf[i],
                     seq.ψ[i],
                     seq.ADC[i],
                     seq.excitation_bool[intervals],
                     seq.t[i],
                     seq.Δt[intervals])
end
Base.view(seq::DiscreteSequence, i::UnitRange) = @views begin
    intervals = i.start:i.stop-1
    DiscreteSequence(seq.Gx[i],
                     seq.Gy[i],
                     seq.Gz[i],
                     seq.B1[i],
                     seq.Δf[i],
                     seq.ψ[i],
                     seq.ADC[i],
                     seq.excitation_bool[intervals],
                     seq.t[i],
                     seq.Δt[intervals])
end
Base.iterate(seq::DiscreteSequence) = (seq[1], 2)
Base.iterate(seq::DiscreteSequence, i) = (i <= length(seq)) ? (seq[i], i+1) : nothing

# -- 0.4. Event-state predicates --------------------------------------------
is_GR_on(seq::DiscreteSequence) =  any(!iszero, seq.Gx) || any(!iszero, seq.Gy) || any(!iszero, seq.Gz)
is_RF_on(seq::DiscreteSequence) =  any(seq.excitation_bool)
is_ADC_on(seq::DiscreteSequence) = any(seq.ADC)
is_GR_off(seq::DiscreteSequence) =  !is_GR_on(seq)
is_RF_off(seq::DiscreteSequence) =  !is_RF_on(seq)
is_ADC_off(seq::DiscreteSequence) = !is_ADC_on(seq)

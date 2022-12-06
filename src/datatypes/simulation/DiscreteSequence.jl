struct DiscreteSequence{T<:Real}
    Gx::AbstractVector{T}
    Gy::AbstractVector{T}
    Gz::AbstractVector{T}
    B1::AbstractVector{Complex{T}}
    Δf::AbstractVector{T}
    ADC::AbstractVector{Bool}
    t::AbstractVector{T}
    Δt::AbstractVector{T}
end

DiscreteSequence(xy::T, z::T) where {T<:Real} = Mag([complex(xy)], [z])
Base.length(seq::DiscreteSequence) = length(seq.Δt)
Base.getindex(seq::DiscreteSequence, i::Integer) = begin
    DiscreteSequence(seq.Gx[i, :],
                     seq.Gy[i, :],
                     seq.Gz[i, :],
                     seq.B1[i, :],
                     seq.Δf[i, :],
                     seq.ADC[i, :],
                     seq.t[i, :],
                     seq.Δt[i, :])
end
Base.getindex(seq::DiscreteSequence, i::UnitRange) = begin
    DiscreteSequence(seq.Gx[i.start:i.stop+1],
                     seq.Gy[i.start:i.stop+1],
                     seq.Gz[i.start:i.stop+1],
                     seq.B1[i.start:i.stop+1],
                     seq.Δf[i.start:i.stop+1],
                     seq.ADC[i],
                     seq.t[i.start:i.stop+1],
                     seq.Δt[i])
end
Base.view(seq::DiscreteSequence, i::UnitRange) = begin
    @views DiscreteSequence(seq.Gx[i.start:i.stop+1],
                     seq.Gy[i.start:i.stop+1],
                     seq.Gz[i.start:i.stop+1],
                     seq.B1[i.start:i.stop+1],
                     seq.Δf[i.start:i.stop+1],
                     seq.ADC[i],
                     seq.t[i.start:i.stop+1],
                     seq.Δt[i])
end
Base.iterate(seq::DiscreteSequence) = (seq[1], 2)
Base.iterate(seq::DiscreteSequence, i) = (i <= length(seq)) ? (seq[i], i+1) : nothing

is_GR_on(seq::DiscreteSequence) =  sum(abs.([seq.Gx; seq.Gy; seq.Gz])) != 0
is_RF_on(seq::DiscreteSequence) =  sum(abs.(seq.B1)) != 0
is_ADC_on(seq::DiscreteSequence) = sum(abs.(seq.ADC)) != 0
is_GR_off(seq::DiscreteSequence) =  !is_GR_on(seq)
is_RF_off(seq::DiscreteSequence) =  !is_RF_on(seq)
is_ADC_off(seq::DiscreteSequence) = !is_ADC_on(seq)
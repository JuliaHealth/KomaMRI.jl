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

Base.length(seq::DiscreteSequence) = length(seq.t)
Base.getindex(seq::DiscreteSequence, i) = DiscreteSequence(
                                                        seq.Gx[i],
                                                        seq.Gy[i],
                                                        seq.Gz[i],
                                                        seq.B1[i],
                                                        seq.Δf[i],
                                                        seq.ADC[i],
                                                        seq.t[i],
                                                        seq.Δt[i]
                                                        )
# Base.view(seq::DiscreteSequence, i) = @views seq[i]
Base.iterate(seq::DiscreteSequence) = (seq[1], 2)
Base.iterate(seq::DiscreteSequence, i) = (i <= length(seq)) ? (seq[i], i+1) : nothing

is_GR_on(seq::DiscreteSequence) =  sum(abs.([seq.Gx; seq.Gy; seq.Gz])) != 0
is_RF_on(seq::DiscreteSequence) =  sum(abs.(seq.B1)) != 0
is_ADC_on(seq::DiscreteSequence) = sum(abs.(seq.ADC)) != 0
is_GR_off(seq::DiscreteSequence) =  !is_GR_on(seq)
is_RF_off(seq::DiscreteSequence) =  !is_RF_on(seq)
is_ADC_off(seq::DiscreteSequence) = !is_ADC_on(seq)
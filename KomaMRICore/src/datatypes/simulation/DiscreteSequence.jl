struct DiscreteSequence{T<:Real}
    Gx::AbstractVector{T}
    Gy::AbstractVector{T}
    Gz::AbstractVector{T}
    B1::AbstractVector{Complex{Float32}}
    Δf::AbstractVector{Complex{Float32}}
    ADC::AbstractVector{Bool}
    t::AbstractVector{T}
    Δt::AbstractVector{T}
end

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

is_GR_on(seq::DiscreteSequence) =  sum(abs.([seq.Gx[1:end-1]; seq.Gy[1:end-1]; seq.Gz[1:end-1]])) != 0
is_RF_on(seq::DiscreteSequence) =  sum(abs.(seq.B1[1:end-1])) != 0
is_ADC_on(seq::DiscreteSequence) = sum(abs.(seq.ADC[1:end-1])) != 0
is_GR_off(seq::DiscreteSequence) =  !is_GR_on(seq)
is_RF_off(seq::DiscreteSequence) =  !is_RF_on(seq)
is_ADC_off(seq::DiscreteSequence) = !is_ADC_on(seq)

function discretize(seq::Sequence; simParams=default_sim_params())
    t, Δt      = get_uniform_times(seq, simParams["Δt"]; Δt_rf=simParams["Δt_rf"])
    B1, Δf     = get_rfs(seq, t)
    Gx, Gy, Gz = get_grads(seq, t)
    tadc       = get_adc_sampling_times(seq)
    ADCflag    = [any(tt .== tadc) for tt in t] #Displaced 1 dt, sig[i]=S(ti+dt)
    seqd       = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    return seqd
end

# Returns the phase in radians
function phase(t, Δt, Δf)
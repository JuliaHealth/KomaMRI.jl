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

"""
Return the DiscreteSequence object with the important vectors prepared for simulation
"""
function discretize(seq::Sequence; simParams=default_sim_params(), isnew=false)
    # To BE DEPRECATED?
    ### Previous discretization
    ### The new discretization is faster but is has more allocation
    ### When simulating, the new simulation with the new discretization is slower, I don't know why ...
    ### ... using ProfileView and Cthulhu can be really helpful for debbuging these perferomances issues
    if !isnew
        t, Δt      = get_uniform_times(seq, simParams["Δt"]; Δt_rf=simParams["Δt_rf"])
        B1, Δf     = get_rfs(seq, t)
        Gx, Gy, Gz = get_grads(seq, t)
        tadc       = get_adc_sampling_times(seq)
        ADCflag    = [any(tt .== tadc) for tt in t]  #Displaced 1 dt, sig[i]=S(ti+dt)
        return DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    end
    # This is the new discretization
    Δtgr, Δtrf = simParams["Δt"], simParams["Δt_rf"]
    sq = samples(seq, Δtgr, Δtrf)
    t, Δt, B1, Δf, Gx, Gy, Gz, ADCflag = sq.t, sq.Δt, sq.rfa, sq.rfΔfc, sq.gxa, sq.gya, sq.gza, sq.adconmask
    return DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
end


############################################################################################
### TEST DISCRETE SEQUENCE STRUCT ##########################################################
############################################################################################
struct SEQD{T<:Real}
    Δt::AbstractVector{T}
    t::AbstractVector{T}
    rfa::AbstractVector{Complex{T}}
    rfΔfc::AbstractVector{T}
    gxa::AbstractVector{T}
    gya::AbstractVector{T}
    gza::AbstractVector{T}
    adconmask::AbstractVector{Bool}
end

function Base.length(seqd::SEQD)
    return length(seqd.Δt)
end
function Base.getindex(seqd::SEQD, i::Integer)
    return SEQD(seqd.Δt[i, :], seqd.t[i, :], seqd.rfa[i, :], seqd.rfΔfc[i, :], seqd.gxa[i, :], seqd.gya[i, :], seqd.gza[i, :], seqd.adconmask[i, :])
end
function Base.getindex(seqd::SEQD, i::UnitRange)
    r = (i.start:i.stop+1)
    return SEQD(seqd.Δt[i], seqd.t[r], seqd.rfa[r], seqd.rfΔfc[r], seqd.gxa[r], seqd.gya[r], seqd.gza[r], seqd.adconmask[r])
end
function Base.view(seqd::SEQD, i::UnitRange)
    r = (i.start:i.stop+1)
    return @views SEQD(seqd.Δt[i], seqd.t[r], seqd.rfa[r], seqd.rfΔfc[r], seqd.gxa[r], seqd.gya[r], seqd.gza[r], seqd.adconmask[r])
end

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

############################################################################################
### For Block and Sequence #################################################################
############################################################################################
"""
Get samples of a block-sequence
"""
function blockvalues(seq::Sequence, blk::Int64; Δtgr::Float64=1e-3, Δtrf::Float64=1e-5)

    # Select the block of the sequence and the events
    #s =
    ΔT = KomaMRICore.durs(seq)[blk]
    rf, gx, gy, gz, adc = seq.RF[blk], seq.GR[1,blk], seq.GR[2,blk], seq.GR[3,blk], seq.ADC[blk]

    # Get the critical times that define the events
    rf_is0, rf_tu, _, _ = KomaMRICore.eventvalues(rf)
    gx_is0, gx_tu, _ = KomaMRICore.eventvalues(gx)
    gy_is0, gy_tu, _ = KomaMRICore.eventvalues(gy)
    gz_is0, gz_tu, _ = KomaMRICore.eventvalues(gz)
    adc_is0, adc_t = KomaMRICore.eventvalues(adc)

    # Get the simulations times
    # This can be optimized.
    # So far, if the Δtrf, Δtgr and Δtadc (an ADC with uniform times, Δtadc is not defined till now)
    # all have a minimum common multiplier equal to the minimum([Δtrf, Δtgr, Δtadc]) = Δtmin,
    # (i.e. (Δtrf = a * Δtmin) and (Δtgr = b * Δtmin) and (Δtadc = c * Δtmin) where any of {a,b,c} is equal to 1)
    # then the simulator doesn't generate unnecessary simulation times

    # Get the uniform-times for RFs and GRs
    toffset = (adc_is0) ? (0.) : (adc_t[1])                      # align with the first sample of the ADC
    trf = [reverse(toffset-Δtrf:-Δtrf:0.); (toffset:Δtrf:ΔT)]    # uniform-times for RFs (they are gonna be included in the simulation when rf is on)
    tgr = [reverse(toffset-Δtgr:-Δtgr:0.); (toffset:Δtgr:ΔT)]    # uniform-times for GRs (they are gonna be included in the simulation when gr is on)

    # Get the simulation-times
    t = Float64[]
    (!rf_is0) && append!(t, trf[(rf_tu[1] .< trf) .& (trf .< rf_tu[end])], [rf_tu[1]; KomaMRICore.get_RF_center(rf); rf_tu[end]])  # consider RF pivot times and equispaced RF times when is RF on
    (!gx_is0) && append!(t, tgr[(gx_tu[1] .< tgr) .& (tgr .< gx_tu[end])], [gx_tu[1]; gx_tu[2]; gx_tu[end-1]; gx_tu[end]])         # consider GX pivot times and equispaced GR times when is GX on
    (!gy_is0) && append!(t, tgr[(gy_tu[1] .< tgr) .& (tgr .< gy_tu[end])], [gy_tu[1]; gy_tu[2]; gy_tu[end-1]; gy_tu[end]])         # consider GY pivot times and equispaced GR times when is GY on
    (!gz_is0) && append!(t, tgr[(gz_tu[1] .< tgr) .& (tgr .< gz_tu[end])], [gz_tu[1]; gz_tu[2]; gz_tu[end-1]; gz_tu[end]])         # consider GZ pivot times and equispaced GR times when is GZ on
    (!adc_is0) && append!(t, adc_t)                                                                  # consider ADC sampling-times
    sort!(t); unique!(t)        # make the simulation-times unique and in increasing order

    # Return for no times
    (isempty(t)) && return [], NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN

    # Get the sampled values of RF and GRs at times t
    rf_a, rf_Δf = KomaMRICore.eventvalues(rf, t)
    gx_a, gy_a, gz_a = KomaMRICore.eventvalues(gx, t), KomaMRICore.eventvalues(gy, t), KomaMRICore.eventvalues(gz, t)
    amps = [rf_a rf_Δf gx_a gy_a gz_a]

    # Create a matrix to add additional amplitudes and times when there are more samples at a certain time
    m = [rf_a rf_Δf gx_a gy_a gz_a]
    tc = Float64[]
    for i in eachindex(t)
        Nsamp = maximum(length.(m[i,:]))
        append!(tc, t[i].*ones(Nsamp))
        if Nsamp > 1
            for j in 1:5
                while length(m[i,j]) < Nsamp
                    push!(m[i,j], m[i,j][end])
                end
            end
        end
    end
    Δtc = tc[2:end] - tc[1:end-1]
    #adc_flag = [any(ti .== adc_t) for ti in tc]
    adc_flag = [any(ti .== adc_t[2:end-1]) for ti in tc]
    if !adc_is0
        adc_flag[findlast(adc_t[1] .== tc)] = true
        adc_flag[findfirst(adc_t[end] .== tc)] = true
    end

    # Return at amplitudes and times with possibly more samples at the same time
    rfa, rfΔf, gxa, gya, gza = vcat(m[:,1]...), vcat(m[:,2]...), vcat(m[:,3]...), vcat(m[:,4]...), vcat(m[:,5]...)
    return tc, Δtc, rfa, rfΔf, gxa, gya, gza, adc_flag, adc_t

end

"""
Get samples of the complete sequence
"""
function sequencevalues(seq::Sequence; Δtgr::Float64=1e-3, Δtrf::Float64=1e-5)
    t, Δt, rfa, rfΔf, gxa, gya, gza, adcflag, adct = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Bool[], Float64[]
    ΔT = durs(seq)
    T0 = cumsum([0; ΔT[:]]) #Start time of each block
    for i = 1:length(seq)
        ti, Δti, rfai, rfΔfi, gxai, gyai, gzai, adcflagi, adcti = blockvalues(seq, i; Δtgr, Δtrf)
        append!(t, T0[i] .+ ti); append!(Δt, Δti)
        append!(rfa, rfai); append!(rfΔf, rfΔfi)
        append!(gxa, gxai); append!(gya, gyai); append!(gza, gzai)
        append!(adcflag, adcflagi); append!(adct, T0[i] .+ adcti)
    end
    return t, Δt, rfa, rfΔf, gxa, gya, gza, adcflag, adct
end

"""
Same as discretize, but returns all the important vector values (not the DiscreteSequence object)
This is for testing purpose only
"""
function seqvals(seq::Sequence; simParams=default_sim_params(), isnew=false)
    if isnew
        t, Δt, B1, Δf, Gx, Gy, Gz, ADCflag, tadc = sequencevalues(seq; Δtgr=simParams["Δt"], Δtrf=simParams["Δt_rf"])
        #ADCflag = [false for _ in eachindex(t)]
        #[ADCflag[i]=true for i in eachindex(tadc)]
        return t, Δt, B1, Δf, Gx, Gy, Gz, ADCflag, tadc
    end
    t, Δt      = get_uniform_times(seq, simParams["Δt"]; Δt_rf=simParams["Δt_rf"])
    B1, Δf     = get_rfs(seq, t)
    Gx, Gy, Gz = get_grads(seq, t)
    tadc       = get_adc_sampling_times(seq)
    ADCflag    = [any(tt .== tadc) for tt in t]  #Displaced 1 dt, sig[i]=S(ti+dt)
    #ADCflag = [false for _ in eachindex(t)]
    #[ADCflag[i]=true for i in eachindex(tadc)]
    return t, Δt, B1, Δf, Gx, Gy, Gz, ADCflag, tadc
end

"""
Return the DiscreteSequence object with the important vectors prepared for simulation
"""
function discretize(seq::Sequence; simParams=default_sim_params(), isnew=false)
    t, Δt, B1, Δf, Gx, Gy, Gz, ADCflag, tadc = seqvals(seq; simParams, isnew)
    return DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
end

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
function discretize(seq::Sequence; simParams=default_sim_params())
    t, Δt      = get_uniform_times(seq, simParams["Δt"]; Δt_rf=simParams["Δt_rf"])
    B1, Δf     = get_rfs(seq, t)
    Gx, Gy, Gz = get_grads(seq, t)
    tadc       = get_adc_sampling_times(seq)
    ADCflag    = [any(tt .== tadc) for tt in t]  #Displaced 1 dt, sig[i]=S(ti+dt)
    return DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
end

############################################################################################
### For Block and Sequence #################################################################
############################################################################################
"""
Get samples of a block-sequence
"""
function blockvalues(seq::Sequence, blk::Int64; Δtgr::Float64=1e-3, Δtrf::Float64=1e-5)

    # Select the block of the sequence and the events
    #s =
    ΔT = durs(seq)[blk]
    rf, gx, gy, gz, adc = seq.RF[blk], seq.GR[1,blk], seq.GR[2,blk], seq.GR[3,blk], seq.ADC[blk]

    # Get the critical times that define the events
    rf_is0, rf_tu, _, _ = eventvalues(rf)
    gx_is0, gx_tu, _ = eventvalues(gx)
    gy_is0, gy_tu, _ = eventvalues(gy)
    gz_is0, gz_tu, _ = eventvalues(gz)
    adc_is0, adc_t = eventvalues(adc)

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
    (!rf_is0) && append!(t, trf[(rf_tu[1] .< trf) .& (trf .< rf_tu[end])], [rf_tu[1]; get_RF_center(rf); rf_tu[end]])  # consider RF pivot times and equispaced RF times when is RF on
    (!gx_is0) && append!(t, tgr[(gx_tu[1] .< tgr) .& (tgr .< gx_tu[end])], [gx_tu[1]; gx_tu[2]; gx_tu[end-1]; gx_tu[end]])         # consider GX pivot times and equispaced GR times when is GX on
    (!gy_is0) && append!(t, tgr[(gy_tu[1] .< tgr) .& (tgr .< gy_tu[end])], [gy_tu[1]; gy_tu[2]; gy_tu[end-1]; gy_tu[end]])         # consider GY pivot times and equispaced GR times when is GY on
    (!gz_is0) && append!(t, tgr[(gz_tu[1] .< tgr) .& (tgr .< gz_tu[end])], [gz_tu[1]; gz_tu[2]; gz_tu[end-1]; gz_tu[end]])         # consider GZ pivot times and equispaced GR times when is GZ on
    (!adc_is0) && append!(t, adc_t)                                                                  # consider ADC sampling-times
    push!(t, ΔT)                # add the last time point always
    sort!(t); unique!(t)        # make the simulation-times unique and in increasing order

    # Return for no times
    (isempty(t)) && return [], NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN

    # Get the sampled values of RF and GRs at times t
    rf_a, rf_Δf = eventvalues(rf, t)
    gx_a, gy_a, gz_a = eventvalues(gx, t), eventvalues(gy, t), eventvalues(gz, t)

    # Create a matrix to add additional amplitudes and times when there are more samples at a certain time
    # so we can get all the times to be simulated (some times could be repeated up to twice)
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

    # Get the masks when the events are on
    Ntc = length(tc)
    rf_onmask = (rf_is0) ? BitVector(zeros(Ntc)) : ((rf_tu[1] .<= tc) .& (tc .<= rf_tu[end]))
    gx_onmask = (gx_is0) ? BitVector(zeros(Ntc)) : ((gx_tu[1] .<= tc) .& (tc .<= gx_tu[end]))
    gy_onmask = (gy_is0) ? BitVector(zeros(Ntc)) : ((gy_tu[1] .<= tc) .& (tc .<= gy_tu[end]))
    gz_onmask = (gz_is0) ? BitVector(zeros(Ntc)) : ((gz_tu[1] .<= tc) .& (tc .<= gz_tu[end]))
    adc_onmask = (adc_is0) ? BitVector(zeros(Ntc)) : ((adc_t[1] .<= tc) .& (tc .<= adc_t[end]))

    # Return at amplitudes and times with possibly more samples at the same time
    rfa, rfΔf, gxa, gya, gza = vcat(m[:,1]...), vcat(m[:,2]...), vcat(m[:,3]...), vcat(m[:,4]...), vcat(m[:,5]...)
    return tc, Δtc, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, adc_t
end

"""
Get samples of the complete sequence
"""
function sequencevalues(seq::Sequence; Δtgr::Float64=1e-3, Δtrf::Float64=1e-5)
    t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, adct = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Bool[], Bool[], Bool[], Bool[], Bool[], Float64[]
    to = cumsum([0; durs(seq)])
    for k = 1:length(seq)
        tk, Δtk, rfak, rfΔfk, gxak, gyak, gzak, rf_onmaskk, gx_onmaskk, gy_onmaskk, gz_onmaskk, adc_onmaskk, adctk = blockvalues(seq, k; Δtgr, Δtrf)
        append!(t, to[k] .+ tk); append!(Δt, Δtk)
        append!(rfa, rfak); append!(rfΔf, rfΔfk)
        append!(gxa, gxak); append!(gya, gyak); append!(gza, gzak)
        append!(rf_onmask, rf_onmaskk); append!(gx_onmask, gx_onmaskk); append!(gy_onmask, gy_onmaskk); append!(gz_onmask, gz_onmaskk); append!(adc_onmask, adc_onmaskk)
        append!(adct, to[k] .+ adctk)
    end
    return t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, adct
end

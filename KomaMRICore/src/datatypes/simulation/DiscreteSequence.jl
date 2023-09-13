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
    ### Previous discretization
    ### The new discretization is faster but is has more allocation
    ### When simulating, the new simulation with the new discretization is slower, I don't know why ...
    ### ... using ProfileView and Cthulhu can be really helpful for debbuging these perferomances issues
    #t, Δt      = get_uniform_times(seq, simParams["Δt"]; Δt_rf=simParams["Δt_rf"])
    #B1, Δf     = get_rfs(seq, t)
    #Gx, Gy, Gz = get_grads(seq, t)
    #tadc       = get_adc_sampling_times(seq)
    #ADCflag    = [any(tt .== tadc) for tt in t]  #Displaced 1 dt, sig[i]=S(ti+dt)
    sq = sequencevalues(seq; Δtgr=simParams["Δt"] , Δtrf=simParams["Δt_rf"])
    t, Δt, B1, Δf, Gx, Gy, Gz, ADCflag = sq.t, sq.Δt, sq.rfa, sq.rfΔf, sq.gxa, sq.gya, sq.gza, sq.adc_onmask
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

    # Get the unique-increasing-times
    t = Float64[]
    (!rf_is0) && append!(t, trf[(rf_tu[1] .< trf) .& (trf .< rf_tu[end])], [rf_tu[1]; get_RF_center(rf); rf_tu[end]])  # consider RF pivot times and equispaced RF times when is RF on
    (!gx_is0) && append!(t, tgr[(gx_tu[1] .< tgr) .& (tgr .< gx_tu[end])], [gx_tu[1]; gx_tu[2]; gx_tu[end-1]; gx_tu[end]])         # consider GX pivot times and equispaced GR times when is GX on
    (!gy_is0) && append!(t, tgr[(gy_tu[1] .< tgr) .& (tgr .< gy_tu[end])], [gy_tu[1]; gy_tu[2]; gy_tu[end-1]; gy_tu[end]])         # consider GY pivot times and equispaced GR times when is GY on
    (!gz_is0) && append!(t, tgr[(gz_tu[1] .< tgr) .& (tgr .< gz_tu[end])], [gz_tu[1]; gz_tu[2]; gz_tu[end-1]; gz_tu[end]])         # consider GZ pivot times and equispaced GR times when is GZ on
    (!adc_is0) && append!(t, adc_t)     # consider ADC sampling-times
    append!(t, [0.; ΔT])                # add the first and last time point always (this is very important when putting blocks together in a sequence)
    sort!(t); unique!(t)                # make the unique-increasing-times actually unique and in increasing order

    # Return for no times
    (isempty(t)) && return Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], BitVector(), BitVector(), BitVector(), BitVector(), BitVector(), Float64[]

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
    adc_onmask = BitVector([any(t .== adc_t) for t in tc])

    # Return amplitudes and simulation-times with possibly more samples at the same time
    rfa, rfΔf, gxa, gya, gza = vcat(m[:,1]...), vcat(m[:,2]...), vcat(m[:,3]...), vcat(m[:,4]...), vcat(m[:,5]...)
    return (t = tc, Δt = Δtc, rfa = rfa, rfΔf = rfΔf, gxa = gxa, gya = gya, gza = gza,
            rf_onmask = rf_onmask, gx_onmask = gx_onmask, gy_onmask = gy_onmask, gz_onmask = gz_onmask,
            adc_onmask = adc_onmask, adc_t = adc_t)
end

"""
Get samples of the complete sequence
"""
function sequencevalues(seq::Sequence; Δtgr::Float64=1e-3, Δtrf::Float64=1e-5)

    # Create empty vectors to be filled
    t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, tadc = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Bool[], Bool[], Bool[], Bool[], Bool[], Float64[]

    # These are the initial times of the blocks
    to = cumsum([0; durs(seq)])

    # Iterate over each block of the sequence
    Nblk = length(seq)                                   # number of blocks
    blk_range = Vector{UnitRange{Int64}}(undef, Nblk)    # empty vector of block-range to be filled
    ko = 1                                               # initialization of the index of the first time sample for block ranges
    for k in 1:Nblk

        # Get the vector values of the block
        blk = blockvalues(seq, k; Δtgr, Δtrf)

        # Fill the vector of block masks
        kn = ko + length(blk.t) - 1     # index of the last time sample of this 1-block-sequence (it is also the first time sample of the next 1-block-sequence)
        blk_range[k] = (ko:kn)          # range of the 1-block-sequence for time samples
        ko = kn                         # update the index of the first time sample for the next block range

        # Add the initial condition just for the first block
        if k == 1
            append!(t, to[k] .+ blk.t[1])
            append!(rfa, blk.rfa[1]); append!(rfΔf, blk.rfΔf[1]); append!(gxa, blk.gxa[1]); append!(gya, blk.gya[1]); append!(gza, blk.gza[1])
            append!(rf_onmask, blk.rf_onmask[1]); append!(gx_onmask, blk.gx_onmask[1]); append!(gy_onmask, blk.gy_onmask[1]); append!(gz_onmask, blk.gz_onmask[1]); append!(adc_onmask, blk.adc_onmask[1])
        end

        # Fill the vector values of the sequence without considering the initial condition,
        # except for the delta-times (the last time point for the sequence duration is always added)
        # and except for the adc-times (wich shouldn't be used, it should be used adc-on-mask instead)
        append!(t, to[k] .+ blk.t[2:end]); append!(Δt, blk.Δt); append!(tadc, to[k] .+ blk.adc_t)
        append!(rfa, blk.rfa[2:end]); append!(rfΔf, blk.rfΔf[2:end]); append!(gxa, blk.gxa[2:end]); append!(gya, blk.gya[2:end]); append!(gza, blk.gza[2:end])
        append!(rf_onmask, blk.rf_onmask[2:end]); append!(gx_onmask, blk.gx_onmask[2:end]); append!(gy_onmask, blk.gy_onmask[2:end]); append!(gz_onmask, blk.gz_onmask[2:end]); append!(adc_onmask, blk.adc_onmask[2:end])

        # Add for the previous-block-last-sample for the masks
        rf_onmask[blk_range[k][1]] |= blk.rf_onmask[1]; gx_onmask[blk_range[k][1]] |= blk.gx_onmask[1]; gy_onmask[blk_range[k][1]] |= blk.gy_onmask[1]; gz_onmask[blk_range[k][1]] |= blk.gz_onmask[1]; adc_onmask[blk_range[k][1]] |= blk.adc_onmask[1]

    end

    # Return the vector values of the sequence
    return (t = t, Δt = Δt, rfa = rfa, rfΔf = rfΔf, gxa = gxa, gya = gya, gza = gza,
            rf_onmask = rf_onmask, gx_onmask = gx_onmask, gy_onmask = gy_onmask, gz_onmask = gz_onmask,
            adc_onmask = adc_onmask, tadc = tadc, blk_range = blk_range)
end

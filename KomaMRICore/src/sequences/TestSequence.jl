using ..KomaMRICore

function check_sequence(seq::Sequence, sys::Scanner; dt=1e-3)

    # Create empty arrays to fill data ...
    Nblk = length(seq)

    # ... for RF hardware constraints
    rf_max   = sys.B1 * ones(Nblk)
    rf_delay = sys.RF_dead_time_T * ones(Nblk)
    rf_Δt    = sys.RF_ring_down_T * ones(Nblk)
    s_rf_past, T0f_rf_past = Sequence([Grad(0, 0)]), 0

    # ... for GR hardware constraints
    gi_max  = sys.Gmax * ones(Nblk, 3)
    gr_max  = sys.Gmax * ones(Nblk)
    sri_max = sys.Smax * ones(Nblk, 3)
    srr_max = sys.Smax * ones(Nblk)

    # ... for ADC hardware constraints
    adc_delay = sys.ADC_dead_time_T * ones(Nblk)

    # Iterate over the blocks of the sequence
    T0 = cumsum([0; durs(seq)])
    for i = 1:Nblk
        T0i, T0f  =  T0[i], T0[i+1]
        s = seq[i]
        # When the RF is on in the block ...
        if is_RF_on(s)
            # ... fill the maximum values of the RF amplitudes
            rf_max[i] = maximum(abs.(s.RF[1].A))
            # ... fill the delays of each RF
            rf_delay[i] = s.RF.delay[1]
            # ... fill the time between two consecutive RFs
            Δt_down_rf_past = s_rf_past[1].DUR[1] - s_rf_past.RF.dur[1]
            Δt_rf_seqs = T0i - T0f_rf_past
            s_rf_past, T0f_rf_past = s, T0f
            rf_Δt[i] = Δt_down_rf_past + Δt_rf_seqs + s.RF.delay[1]
        end
        # When any of the GR is on in the block ...
        if is_GR_on(s)
            # ... fill the maximum amplitudes of every GR x,y,z component
            gi_max[i, :] = [maximum(abs.(gi.A)) for gi in s.GR]
            # ... fill the maximum module values for each sample of the GR
            t, Δt = get_uniform_times(s, dt)
            gx, gy, gz = get_grads(s, t)
            gr_max[i] = maximum(sqrt.(gx.^2 + gy.^2 + gz.^2))
            # ... fill the maximum slew-rates of every GR x,y,z component
            sri_ini = [gi.rise == 0 && gi.A[1]   == 0 ? 0 : abs(gi.A[1]   / gi.rise) for gi in s.GR]
            sri_end = [gi.fall == 0 && gi.A[end] == 0 ? 0 : abs(gi.A[end] / gi.fall) for gi in s.GR]
            sri_mid = [length(gi.A) == 1 ? 0 : maximum(abs.((gi.A[2:end] - gi.A[1:end-1]) ./ (length(gi.T) == 1 ? gi.T / (length(gi.A) - 1) : gi.T))) for gi in s.GR]
            sri_max[i, :] = maximum([sri_ini sri_mid sri_end], dims=2)
            # ... fill the maximum module values for each sample of the GR slew-rates
            srx, sry, srz = (gx[2:end] - gx[1:end-1]) ./ Δt, (gy[2:end] - gy[1:end-1]) ./ Δt, (gz[2:end] - gz[1:end-1]) ./ Δt
            srr_max[i] = maximum(sqrt.(srx.^2 + sry.^2 + srz.^2))
        end
        # When the ADC is on in the block ...
        if is_ADC_on(s)
            # ... fill the delays of each ADC
            adc_delay[i] = s.ADC.delay[1]
        end
    end

    # Generate error tables for magnitudes and print info to user
    mag_err_cols = ["B1"; "GX"; "GY"; "GZ"; "GR"; "SX"; "SY"; "SZ"; "SR"]
    mag_err_all = [rf_max .> sys.B1 gi_max .> sys.Gmax gr_max .> sys.Gmax sri_max .> sys.Smax srr_max .> sys.Smax]
    mag_err_any = [any(mag_err_all, dims=1)...]
    mag_err_blks = findall([any(mag_err_all, dims=2)...])
    mag_err_needed = mag_err_all[mag_err_blks, :]
    mag_error_table = [mag_err_blks mag_err_needed]
    if any(mag_err_any)
        @error "HW magnitude errors in: $(mag_err_cols[findall(mag_err_any)])\n $(["BLK"; mag_err_cols])" mag_error_table
    else
        @info "HW magnitude constraints are satisfied"
    end

    # Generate error tables for timing and print info to user
    time_err_cols = ["RF_dead_time"; "RF_ring_down_time"; "ADC_dead_time"]
    time_err_all = [rf_delay .< sys.RF_dead_time_T rf_Δt .< sys.RF_ring_down_T adc_delay .< sys.ADC_dead_time_T]
    time_err_any = [any(time_err_all, dims=1)...]
    time_err_blks = findall([any(time_err_all, dims=2)...])
    time_err_needed = time_err_all[time_err_blks, :]
    time_error_table = [time_err_blks time_err_needed]
    if any(time_err_any)
        @error "HW timing errors in: $(time_err_cols[findall(time_err_any)])\n $(["BLK"; time_err_cols])" time_error_table
    else
        @info "HW timing constraints are satisfied"
    end

    # Return the error tables
    return mag_error_table, time_error_table
end

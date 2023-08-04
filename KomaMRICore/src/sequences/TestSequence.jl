using ..KomaMRICore

function check_sequence(seq::Sequence, sys::Scanner; dt=1e-3)
    is_hw_compliant = true
    rf_max_amps = Float64[]
    gr_max_amps = reshape(Float64[], 3, 0)
    gr_max_abs_amps = Float64[]
    gr_max_slewrates = reshape(Float64[], 3, 0)
    gr_max_abs_slewrates = Float64[]
    s_adc_past, T0f_adc_past = Sequence([Grad(0, 0)]), 0
    Δt_delay_adcs = Float64[]
    Δt_between_adcs = Float64[]
    s_rf_past, T0f_rf_past = Sequence([Grad(0, 0)]), 0
    Δt_delay_rfs = Float64[]
    Δt_between_rfs = Float64[]
    T0 = cumsum([0; durs(seq)])
    for i = 1:length(seq)
        T0i, T0f  =  T0[i], T0[i+1]
        s = seq[i]
        if is_RF_on(s)
            rf = maximum(abs.(s.RF[1].A))
            append!(rf_max_amps, rf)
            #
            append!(Δt_delay_rfs, s.RF.delay[1])
            # Get the
            Δt_rf_past = s_rf_past[1].DUR[1] - s_rf_past.RF.dur[1]
            append!(Δt_between_rfs, s.RF.delay[1] + T0i - T0f_rf_past - Δt_rf_past)
            s_rf_past, T0f_rf_past = s, T0f
        end
        if is_GR_on(s)
            # Get maximum amplitudes of every x,y,z component
            gr = [maximum(abs.(gi.A)) for gi in s.GR]
            gr_max_amps = hcat(gr_max_amps, gr)
            # Get maximum absolute values for each sample of the grads
            t, Δt = get_uniform_times(s, dt)
            gx, gy, gz = get_grads(s, t)
            gr = maximum(sqrt.(gx.^2 + gy.^2 + gz.^2))
            append!(gr_max_abs_amps, gr)
            # Get maximum slew-rates of every x,y,z component
            sr_ini = [gi.rise == 0 && gi.A[1]   == 0 ? 0 : abs(gi.A[1]   / gi.rise) for gi in s.GR]
            sr_end = [gi.fall == 0 && gi.A[end] == 0 ? 0 : abs(gi.A[end] / gi.fall) for gi in s.GR]
            sr_mid = [length(gi.A) == 1 ? 0 : maximum(abs.((gi.A[2:end] - gi.A[1:end-1]) ./ (length(gi.T) == 1 ? gi.T / (length(gi.A) - 1) : gi.T))) for gi in s.GR]
            gr_max_slewrates = hcat(gr_max_slewrates, sr_ini, sr_mid, sr_end)
            # Get maximum absolute values for each sample of the slew-rates
            srx, sry, srz = (gx[2:end] - gx[1:end-1]) ./ Δt, (gy[2:end] - gy[1:end-1]) ./ Δt, (gz[2:end] - gz[1:end-1]) ./ Δt
            sr = maximum(sqrt.(srx.^2 + sry.^2 + srz.^2))
            append!(gr_max_abs_slewrates, sr)
        end
        if is_ADC_on(s)
            append!(Δt_delay_adcs, s.ADC.delay[1])
            # Get the
            Δt_adc_past = s_adc_past[1].DUR[1] - s_adc_past.ADC.dur[1]
            append!(Δt_between_adcs, s.ADC.delay[1] + T0i - T0f_adc_past - Δt_adc_past)
            s_adc_past, T0f_adc_past = s, T0f
        end
    end
    if maximum(rf_max_amps) > sys.B1
        @warn "A B1 sample of the sequence is greater than the maximum value B1 of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_amps) > sys.Gmax
        @warn "A GR sample of the sequence is greater than the maximum value Gmax of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_abs_amps) > sys.Gmax
        @warn "An oblique GR sample of the sequence is greater than the maximum value Gmax of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_slewrates) > sys.Smax
        @warn "A GR slew-rate of the sequence is greater than the maximum value Smax of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_abs_slewrates) > sys.Smax
        @warn "An oblique GR slew-rate of the sequence is greater than the maximum value Smax of the scanner"
        is_hw_compliant = false
    end
    if minimum(Δt_delay_adcs) < sys.ADC_dead_time_T
        @warn "The delay of an ADC in a sequence block is lower than the ADC_dead_time_T of the scanner"
        is_hw_compliant = false
    end
    if minimum(Δt_delay_rfs) < sys.RF_dead_time_T
        @warn "The delay of an RF in a sequence block is lower than the RF_dead_time_T of the scanner"
        is_hw_compliant = false
    end
    if minimum(Δt_between_rfs) < sys.RF_ring_down_T
        @warn "A distance between 2 RFs is lower than the RF_ring_down_T of the scanner"
        is_hw_compliant = false
    end
    return Δt_between_adcs
end

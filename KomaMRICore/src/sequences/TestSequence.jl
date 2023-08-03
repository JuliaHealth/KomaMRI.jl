using ..KomaMRICore

function check_sequence(seq::Sequence, sys::Scanner; dt=1e-3)
    is_hw_compliant = true
    rf_max_amps = Float64[]
    gr_max_amps = reshape(Float64[], 3, 0)
    gr_max_abs_amps = Float64[]
    gr_max_slewrates = reshape(Float64[], 3, 0)
    gr_max_abs_slewrates = Float64[]
    for i = 1:2#length(seq)
        s = seq[i]
        if is_RF_on(s)
            rf = maximum(abs.(s.RF[1].A))
            append!(rf_max_amps, rf)
        end
        if is_GR_on(s)
            # Get maximum amplitudes of every x,y,z component
            gx, gy, gz = maximum(abs.(s.GR[1,1].A)), maximum(abs.(s.GR[2,1].A)), maximum(abs.(s.GR[3,1].A))
            gr_max_amps = hcat(gr_max_amps, [gx; gy; gz])
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
    end
    if maximum(rf_max_amps) > sys.B1
        @warn "A B1 sample of the sequence is higher than the maximum value B1 of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_amps) > sys.Gmax
        @warn "A GR sample of the sequence is higher than the maximum value Gmax of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_abs_amps) > sys.Gmax
        @warn "An oblique GR sample of the sequence is higher than the maximum value Gmax of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_slewrates) > sys.Smax
        @warn "A GR slew-rate of the sequence is higher than the maximum value Smax of the scanner"
        is_hw_compliant = false
    end
    if maximum(gr_max_abs_slewrates) > sys.Smax
        @warn "An oblique GR slew-rate of the sequence is higher than the maximum value Smax of the scanner"
        is_hw_compliant = false
    end
    return gr_max_slewrates, gr_max_abs_slewrates
end

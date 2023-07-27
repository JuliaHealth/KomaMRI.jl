using ..KomaMRICore

function check_sequence(seq::Sequence, sys::Scanner)
    is_hw_compliance = true
    rf_max_amps = []
    gr_max_amps = reshape([], 3, 0)
    gr_max_abs_amps = []
    for i = 1:length(seq)
        s = seq[i]
        if is_RF_on(s)
            rf = maximum(abs.(s.RF[1].A))
            append!(rf_max_amps, rf)
        end
        if is_GR_on(s)
            gx, gy, gz = maximum(abs.(s.GR[1,1].A)), maximum(abs.(s.GR[2,1].A)), maximum(abs.(s.GR[3,1].A))
            gr_max_amps = hcat(gr_max_amps, [gx; gy; gz])
            append!(gr_max_abs_amps, sqrt(gx^2 + gy^2 + gz^2))
        end
    end
    if maximum(rf_max_amps) > sys.B1
        @warn "A B1 sample of the sequence is higher than the maximum value B1 of the scanner"
        is_hw_compliance = false
    end
    if maximum(gr_max_amps) > sys.Gmax
        @warn "A GR sample of the sequence is higher than the maximum value Gmax of the scanner"
        is_hw_compliance = false
    end
    if maximum(gr_max_abs_amps) > sys.Gmax
        @warn "An oblique GR sample of the sequence is higher than the maximum value Gmax of the scanner"
        is_hw_compliance = false
    end
    return is_hw_compliance
end

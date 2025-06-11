function generate_unit_sinc_waveform(duration, TBP, sys::Scanner; apodization=0.5)
    Δt = sys.RF_Δt  # Raster time [s]
    n_steps = round(Int, duration / Δt)
    t = range(-duration / 2, stop = duration / 2, length = n_steps + 1)
    bw = TBP / duration

    # Hanning or Hamming window
    window = (1 .- apodization) .+ apodization .* cos.(2π .* collect(-n_steps ÷ 2 : n_steps ÷ 2) ./ n_steps)

    wf = sinc.(bw .* t) .* window
    wf .-= wf[1]  # remove DC offset

    return t, wf
end

function scale_rf_waveform(unit_wf, flip_angle_rad, sys::Scanner)
    Δt = sys.RF_Δt
	gamma_rad = 2 * π * γ

    # Integration with trapezoidal rule
    integral = sum((unit_wf[1:end-1] .+ unit_wf[2:end]) ./ 2) * Δt
    unit_flip_angle = gamma_rad * integral 

    return (flip_angle_rad / unit_flip_angle)
end
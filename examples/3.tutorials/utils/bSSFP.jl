function bSSFP_cine(
    FOV::Float64, 
	N_matrix::Int, 
	TR::Float64, 
	flip_angle, 
	RRs,
	N_phases,
	sys::Scanner; 
	G=[0,0,0], 
	Δf=0,
    adc_duration=1e-3,
    N_dummy_cycles=10
)
    function RR(i)
        idx = (i - 1) % length(RRs) + 1
        return RRs[idx]
    end

    seq = Sequence()
    base_seq =  bSSFP(FOV, N_matrix, TR, flip_angle, sys; Δf=Δf, adc_duration=adc_duration)

    n = 1
    for i in 0:(N_matrix - 1) # 1 vps (Views Per Segment)
        line = base_seq[6*i .+ (1:6)]

        if N_dummy_cycles > 0 && i == 0
            dummy_dur = N_dummy_cycles * dur(line)
            rr_sum = RR(n)
            while dummy_dur > rr_sum
                n += 1
                rr_sum += RR(n)
                
            end
            d_pre = rr_sum - dummy_dur
            seq += Delay(d_pre)
        end
        
        for j in 1:(N_phases + N_dummy_cycles)
            l = copy(line)

            l[1].RF[1].A  *= (-1)^(j)         # Sign of the RF pulse is alteranted every TR
            l[4].ADC[1].ϕ  = j%2==0 ? 0 : π   # so, phase of the ADC is consecuently alterned between 0 and π

            if j in (1:N_dummy_cycles)
                l = 0*l
                l[4].ADC[1].N = 0
            end

            seq += l
        end

        phases_and_dummy_dur = (N_phases + N_dummy_cycles)* dur(line)
        n += 1
        rr_sum = RR(n)
        while phases_and_dummy_dur > rr_sum
            n += 1
            rr_sum += RR(n)
           
        end
        d_post = rr_sum - phases_and_dummy_dur
        seq += Delay(d_post)
    end

    return seq
end

function bSSFP(
	FOV, 
	N, 
	TR, 
	flip_angle, 
	sys; 
	Δf = 0, 
	pulse_duration = 3e-3,
	z0 = 0.0,
	slice_thickness = 10e-3,
	TBP = 4,
	adc_duration = 1e-3
)
	# SINC pulse parameters 
	BW = TBP / pulse_duration

	flip_angle_rad = flip_angle * π / 180 
	t_rf, unit_wf = generate_unit_sinc_waveform(pulse_duration, TBP, sys; apodization=0.5)
	amplitude = scale_rf_waveform(unit_wf, flip_angle_rad, sys)

	# Choose slice thickness
	Gss = BW / (γ * slice_thickness)
	f0 = γ * z0 * Gss

	ζ_ss = Gss / sys.Smax
	area_Gss = Gss * (pulse_duration + ζ_ss)
	T_refocus = (area_Gss/(2*Gss)) - ζ_ss

	phase_dur = T_refocus + 2*ζ_ss

	Δk = (1/FOV)
	FOVk = (N-1)*Δk
	Kmax = FOVk/2
	area = Kmax / γ

	ζ_phase = (-phase_dur * sys.Smax + sqrt((phase_dur * sys.Smax)^2 - 4*sys.Smax*area))/(-2*sys.Smax)
	T_phase = phase_dur - 2*ζ_phase
	G_phase = ζ_phase * sys.Smax

	EX = Sequence([Grad(0., pulse_duration, ζ_ss)  Grad(-G_phase,    T_phase,   ζ_phase); 
				   Grad(0., pulse_duration, ζ_ss)  Grad(-G_phase,    T_phase,   ζ_phase); 
				   Grad(Gss, pulse_duration, ζ_ss) Grad(-Gss,        T_refocus, ζ_ss);], 
				  [RF(amplitude .* unit_wf, pulse_duration, f0, ζ_ss) RF(0,0)])

	# Acquisition ----------------------------------------------
	# PHASE
	area_step = Δk / γ
	G_step = area_step / (T_phase + ζ_phase)

	# FE and Readout
	area_ro = FOVk / γ
	adc_duration_min = (N-1)*sys.ADC_Δt
	if adc_duration < adc_duration_min
		@warn "ADC duration is too short. It will be extended to $(adc_duration_min * 1e3) ms."
		T_ro = adc_duration_min
	else
		T_ro = adc_duration
	end
	ζ_ro = (-T_ro * sys.Smax + sqrt((T_ro * sys.Smax)^2 + 4*sys.Smax*area_ro))/(2*sys.Smax)
	G_ro = ζ_ro * sys.Smax

	ro = Sequence([Grad(G_ro,T_ro,ζ_ro); Grad(0,0); Grad(0,0);;])
	ro.ADC[1] = ADC(N, T_ro, ζ_ro)

	bssfp = Sequence()
	for i in 0:(N-1)
		# Excitation and first phase 
		ex = copy(EX)
		ex[end].GR[2].A += i*G_step

		# FE and Readout
		balance = Sequence([ex[end].GR[1]; -ex[end].GR[2]; ex[end].GR[3];;])

		delay_TR = TR - (dur(ex) + dur(ro) + dur(balance))
		bssfp += (ex + Delay(delay_TR/2) + ro + Delay(delay_TR/2) + balance)
	end
	bssfp.DEF = Dict("Nx"=>N,"Ny"=>N,"Nz"=>1,"Name"=>"bssfp"*string(N)*"x"*string(N),"FOV"=>[FOV, FOV, 0])
	return bssfp

	return EX
end

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
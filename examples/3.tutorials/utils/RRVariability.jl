import KomaMRI.PulseDesigner as PD

function bSSFP_cine(FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    G=[0,0,0], Δf=0, adc_duration=1e-3, N_dummy_cycles=10)
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
            @addblock seq += Delay(d_pre)
        end
        
        for j in 1:(N_phases + N_dummy_cycles)
            l = copy(line)

            l[1].RF[1].A  *= (-1)^(j)         # Sign of the RF pulse is alteranted every TR
            l[4].ADC[1].ϕ  = j%2==0 ? 0 : π   # so, phase of the ADC is consecuently alterned between 0 and π

            if j in (1:N_dummy_cycles)
                l = 0*l
                l[4].ADC[1].N = 0
            end

            @addblock seq += l
        end

        phases_and_dummy_dur = (N_phases + N_dummy_cycles)* dur(line)
        n += 1
        rr_sum = RR(n)
        while phases_and_dummy_dur > rr_sum
            n += 1
            rr_sum += RR(n)
           
        end
        d_post = rr_sum - phases_and_dummy_dur
        @addblock seq += Delay(d_post)
    end

    return seq
end

function bSSFP(FOV, N, TR, flip_angle, sys;
    Δf=0, pulse_duration=3e-3, z0=0.0, slice_thickness=10e-3, TBP=4, adc_duration=1e-3)
    Δk = 1 / FOV
    ro_area = (N - 1) * Δk
    rf, gz, gz_reph = PD.make_sinc_pulse(
        flip_angle * π / 180;
        duration=pulse_duration,
        slice_thickness,
        freq_offset=Δf + TBP / pulse_duration * z0 / slice_thickness,
        time_bw_product=TBP,
        apodization=0.5,
        sys,
    )

    readout_time = max(adc_duration, (N - 1) * sys.ADC_Δt)
    readout_time == adc_duration || @warn "ADC duration is too short. It will be extended to $(readout_time * 1e3) ms."
    gx = PD.make_trapezoid(; flat_area=ro_area, flat_time=readout_time, sys)
    adc = PD.make_adc(N; duration=gx.T, delay=gx.rise, sys)

    excitation_time = max(dur(rf), dur(gz))
    pre_time = dur(gz_reph)
    gx_pre = PD.make_trapezoid(; area=-area(gx) * γ / 2, duration=pre_time, sys)

    bssfp = Sequence(sys)
    @addblock for i in 0:(N - 1)
        ky = i * Δk - ro_area / 2
        gy_pre = PD.make_trapezoid(;
            area=ky,
            duration=pre_time,
            rise_time=gx_pre.rise,
            fall_time=gx_pre.fall,
            sys,
        )
        line_time = excitation_time + pre_time + dur(gx) + pre_time
        delay_TR = TR - line_time
        bssfp += (rf, z=gz)
        bssfp += (x=gx_pre, y=gy_pre, z=gz_reph)
        bssfp += Delay(delay_TR / 2)
        bssfp += (adc, x=gx)
        bssfp += Delay(delay_TR / 2)
        bssfp += (x=gx_pre, y=-gy_pre, z=gz_reph)
    end
    bssfp.DEF = Dict("Nx"=>N, "Ny"=>N, "Nz"=>1, "Name"=>"bssfp$(N)x$(N)", "FOV"=>[FOV, FOV, 0])
    return bssfp
end

function plot_cine(frames, fps; Δt=1/fps, filename="cine_recon.gif", width=400, height=400)

	x = 0:size(frames[1])[2]-1
	y = 1:size(frames[1])[1]

	global_min = minimum(reduce(vcat, frames))
    global_max = maximum(reduce(vcat, frames))

	t = 0

	if ndims(frames) == 1
		n_frames = length(frames)
		n_cines = 1
	else
		n_frames, n_cines = size(frames)
	end

	anim = @animate for i in 1:n_frames
		t += Δt
		l = (1, n_cines)
		plots = []
		for j in 1:n_cines
			push!(plots, Plots.plot!(
				Plots.heatmap(
					x,y,frames[i,j]',color=:greys; 
					aspect_ratio=:equal, 
					colorbar=true, 
					clim=(global_min, global_max),
					size=(n_cines * width, height),
				),
				title="t = "*Printf.@sprintf("%.3f", t)*"s", 
				xlims=(minimum(x), maximum(x)), 
				ylims=(minimum(y), maximum(y))
			))
		end
		Plots.plot(plots..., layout=l)
	end

	gif(anim, filename, fps = fps);
end

function reconstruct_cine(raw, seq, N_matrix, N_phases)
	frames = []
	@info "Running reconstruction"
	@time begin
		recParams = Dict{Symbol,Any}(:reco=>"direct")
		Nx = Ny = N_matrix
		recParams[:reconSize] = (Nx, Ny)
		recParams[:densityWeighting] = false
		acqData = AcquisitionData(raw)
		_, ktraj = get_kspace(seq)
		for i in 1:N_phases
			acqAux = deepcopy(acqData)
			range = reduce(vcat,[j*(N_matrix*N_phases).+((i-1)*N_matrix.+(1:N_matrix)) for j in 0:N_matrix-1])
			## Kdata
			acqAux.kdata[1] = reshape(acqAux.kdata[1][range],(N_matrix^2,1))
			## Traj
			acqAux.traj[1].circular = false
			acqAux.traj[1].nodes = transpose(ktraj[:, 1:2])[:,range]
			acqAux.traj[1].nodes = acqAux.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqAux.traj[1].nodes[:]))
			acqAux.traj[1].numProfiles = N_matrix
			acqAux.traj[1].times = acqAux.traj[1].times[range]
			## Reconstruction
			aux = @timed reconstruction(acqAux, recParams)
			image  = reshape(aux.value.data,Nx,Ny,:)
			image_aux = abs.(image[:,:,1])
			push!(frames,image_aux)
		end
	end
	return frames
end

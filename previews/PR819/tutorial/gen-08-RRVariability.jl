using KomaMRI, PlotlyBase, Printf, Suppressor; #hide
import KomaMRI.PulseDesigner as PD #hide
sys = Scanner(); #hide
function bSSFP_cine(FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; #hide
    Δf=0, adc_duration=1e-3) #hide
    RR(i) = RRs[(i - 1) % length(RRs) + 1] #hide
    base_seq = bSSFP(FOV, N_matrix, TR, flip_angle, sys; Δf, adc_duration) #hide
    seq = Sequence() #hide
    cycle = 0 #hide
    preparation = base_seq[1:6] #hide
    for _ in 1:round(Int, RR(1) / TR) #hide
        cycle += 1 #hide
        line = copy(preparation) #hide
        line[1].RF[1].A *= (-1)^cycle #hide
        line[4].ADC[1].ϕ = iseven(cycle) ? 0 : π #hide
        line[4].ADC[1].N = 0 #hide
        @addblock seq += line #hide
    end #hide
    for i in 0:(N_matrix - 1) #hide
        rr = RR(i + 2) #hide
        line = base_seq[6 * i .+ (1:6)] #hide
        n_cycles = round(Int, rr / TR) #hide
        phase_stride = n_cycles ÷ N_phases #hide
        acquisition_cycles = phase_stride:phase_stride:n_cycles #hide
        for repetition in 1:n_cycles #hide
            cycle += 1 #hide
            l = copy(line) #hide
            l[1].RF[1].A *= (-1)^cycle #hide
            l[4].ADC[1].ϕ = iseven(cycle) ? 0 : π #hide
            repetition in acquisition_cycles || (l[4].ADC[1].N = 0) #hide
            @addblock seq += l #hide
        end #hide
    end #hide
    return seq #hide
end #hide
function bSSFP(FOV, N, TR, flip_angle, sys; #hide
    Δf=0, pulse_duration=3e-3, z0=0.0, slice_thickness=10e-3, TBP=4, adc_duration=1e-3) #hide
    Δk = 1 / FOV #hide
    ro_area = (N - 1) * Δk #hide
    rf, gz, gz_reph = PD.make_sinc_pulse( #hide
        flip_angle * π / 180; #hide
        duration=pulse_duration, #hide
        slice_thickness, #hide
        freq_offset=Δf + TBP / pulse_duration * z0 / slice_thickness, #hide
        time_bw_product=TBP, #hide
        apodization=0.5, #hide
        sys, #hide
    ) #hide
    readout_time = max(adc_duration, (N - 1) * sys.ADC_Δt) #hide
    readout_time == adc_duration || @warn "ADC duration is too short. It will be extended to $(readout_time * 1e3) ms." #hide
    gx = PD.make_trapezoid(; flat_area=ro_area, flat_time=readout_time, sys) #hide
    adc = PD.make_adc(N; duration=gx.T, delay=gx.rise, sys) #hide
    excitation_time = max(dur(rf), dur(gz)) #hide
    pre_time = dur(gz_reph) #hide
    gx_pre = PD.make_trapezoid(; area=-area(gx) * γ / 2, duration=pre_time, sys) #hide
    bssfp = Sequence(sys) #hide
    @addblock for i in 0:(N - 1) #hide
        ky = i * Δk - ro_area / 2 #hide
        gy_pre = PD.make_trapezoid(; area=ky, duration=pre_time, rise_time=gx_pre.rise, fall_time=gx_pre.fall, sys) #hide
        delay_TR = TR - excitation_time - pre_time - dur(gx) - pre_time #hide
        bssfp += (rf, z=gz) #hide
        bssfp += (x=gx_pre, y=gy_pre, z=gz_reph) #hide
        bssfp += Delay(delay_TR / 2) #hide
        bssfp += (adc, x=gx) #hide
        bssfp += Delay(delay_TR / 2) #hide
        bssfp += (x=gx_pre, y=-gy_pre, z=gz_reph) #hide
    end #hide
    bssfp.DEF = Dict("Nx"=>N, "Ny"=>N, "Nz"=>1, "Name"=>"bssfp$(N)x$(N)", "FOV"=>[FOV, FOV, 0]) #hide
    return bssfp #hide
end #hide
function plot_cine(frames, fps; Δt=1 / fps, height=400) #hide
    n_frames, n_cines = ndims(frames) == 1 ? (length(frames), 1) : size(frames) #hide
    zmin, zmax = extrema(reduce(vcat, frames)) #hide
    axis(j) = j == 1 ? "" : string(j) #hide
    trace(i, j) = heatmap(; #hide
        z=permutedims(frames[i, j]), transpose=false, #hide
        xaxis="x$(axis(j))", yaxis="y$(axis(j))", #hide
        coloraxis="coloraxis", #hide
    ) #hide
    frames_plot = [ #hide
        PlotlyBase.frame(; #hide
            name=string(i), data=[trace(i, j) for j in 1:n_cines], #hide
            layout=attr(title=Printf.@sprintf("t = %.3f s", i * Δt)), #hide
        ) #hide
        for i in 1:n_frames #hide
    ] #hide
    animation = attr(; #hide
        mode="immediate", fromcurrent=true, transition=attr(duration=0), #hide
        frame=attr(duration=round(Int, 1000 / fps), redraw=true), #hide
    ) #hide
    layout = Layout(; #hide
        title=Printf.@sprintf("t = %.3f s", Δt), height, #hide
        meta=attr(koma=attr(loop=true)), #hide
        margin=attr(l=40, r=80, t=70, b=40), #hide
        grid=attr(rows=1, columns=n_cines, pattern="independent"), #hide
        coloraxis=attr(cmin=zmin, cmax=zmax, colorscale=[[0, "black"], [1, "white"]]), #hide
        updatemenus=[attr(; #hide
            type="buttons", direction="right", x=0, y=1.1, showactive=false, #hide
            buttons=[ #hide
                attr(label="Play", method="animate", args=[nothing, animation]), #hide
                attr(label="Pause", method="animate", args=[[nothing], attr(mode="immediate", frame=attr(duration=0, redraw=false))]), #hide
            ], #hide
        )], #hide
    ) #hide
    for j in 1:n_cines #hide
        suffix = axis(j) #hide
        layout[Symbol("xaxis$suffix")] = attr(constrain="domain") #hide
        layout[Symbol("yaxis$suffix")] = attr(scaleanchor="x$suffix") #hide
    end #hide
    return Plot( #hide
        [trace(1, j) for j in 1:n_cines], layout, frames_plot; #hide
        config=PlotConfig(displaylogo=false, responsive=true), #hide
    ) #hide
end #hide
function reconstruct_cine(raw, seq, N_matrix, N_phases) #hide
    frames = [] #hide
    recParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(N_matrix, N_matrix), :densityWeighting=>false) #hide
    acqData = AcquisitionData(raw) #hide
    _, ktraj = get_kspace(seq) #hide
    for i in 1:N_phases #hide
        acqAux = deepcopy(acqData) #hide
        range = reduce(vcat, [j * (N_matrix * N_phases) .+ ((i - 1) * N_matrix .+ (1:N_matrix)) for j in 0:(N_matrix - 1)]) #hide
        acqAux.kdata[1] = reshape(acqAux.kdata[1][range], (N_matrix^2, 1)) #hide
        acqAux.traj[1].circular = false #hide
        acqAux.traj[1].nodes = transpose(ktraj[:, 1:2])[:, range] #hide
        acqAux.traj[1].nodes = acqAux.traj[1].nodes[1:2, :] ./ maximum(2 * abs.(acqAux.traj[1].nodes[:])) #hide
        acqAux.traj[1].numProfiles = N_matrix #hide
        acqAux.traj[1].times = acqAux.traj[1].times[range] #hide
        push!(frames, abs.(reshape(reconstruction(acqAux, recParams).data, N_matrix, N_matrix, :)[:, :, 1])) #hide
    end #hide
    return frames #hide
end #hide
nothing #hide

obj = heart_phantom(; rotation_angle=0.0, spins_per_voxel=20);

p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21); #hide
display(p1);

obj.motion

RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 32          # image size = N x N
N_phases     = 16          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 6.25e-3     # [s] bSSFP repetition time
frame_interval = RRs[1] / N_phases # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3;     # [s]

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    adc_duration = adc_duration,
);

# Simulation  #hide
raw1 = @suppress simulate(obj, seq, sys); #hide
# Reconstruction #hide
frames1 = @suppress reconstruct_cine(raw1, seq, N_matrix, N_phases); #hide

fps = 1 / frame_interval #hide
p2 = plot_cine(frames1, fps; Δt=frame_interval); #hide
display(p2);

RRs = [900, 1100, 1000, 1000, 1000, 800] .* 1e-3;

# Take the time curve from the contraction motion:
t_curve = obj.motion.motions[1].time
# Generate a new time curve:
t_curve_new = TimeCurve(
    t = t_curve.t,
    t_unit = t_curve.t_unit,
    periodic = true,
    periods = RRs
)
# Assign the new time curve to both stored motion components:
obj.motion.motions[1].time = t_curve_new
obj.motion.motions[2].time = t_curve_new;

p3 = plot_phantom_map(obj, :T1 ; height=450, time_samples=41); #hide

display(p3);

# Simulation  #hide
raw2 = @suppress simulate(obj, seq, sys) #hide
# Reconstruction #hide
frames2 = @suppress reconstruct_cine(raw2, seq, N_matrix, N_phases); #hide
display(plot_cine(frames2, fps; Δt=frame_interval));

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    adc_duration = adc_duration,
);

# Simulation  #hide
raw3 = @suppress simulate(obj, seq, sys) #hide
# Reconstruction #hide
frames3 = @suppress reconstruct_cine(raw3, seq, N_matrix, N_phases); #hide

display(plot_cine(frames3, fps; Δt=frame_interval));

p4 = plot_cine([frames2 ;; frames3], fps; Δt=frame_interval); #hide
display(p4);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

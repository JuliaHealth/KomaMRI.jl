# # Cardiac Cine MRI with Arrhythmias
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
# This tutorial shows how to simulate cardiac cine MRI using Koma,  
# including cases with variable RR intervals (i.e., arrhythmias). You'll learn how to:

# 1. Simulate a clean cine acquisition with constant RR intervals.  
# 2. Introduce arrhythmias (variable RR intervals) into the cardiac phantom.  
# 3. Observe how this desynchronization degrades image quality.  
# 4. Correct the acquisition by synchronizing the sequence with the phantom’s RR variability (manual triggering).

# ### 1. Constant RR for Phantom and Sequence
#
# We will begin by simulating a cardiac cine on a myocardial phantom with a constant RR interval.
# We'll use the `heart_phantom` function to create a ring-shaped phantom filled with blood,
# resembling the left ventricle. We disable its rotation so this example isolates cardiac
# synchronization from intra-readout rotational motion:

obj = heart_phantom(; rotation_angle=0.0, spins_per_voxel=20);

# The phantom now exhibits periodic contraction with a 1-second period:
p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21); #hide
#md p1 #hide
#jl display(p1);
#nb display(p1);

# As shown in previous tutorials, the phantom's motion is defined by its `motion` field.
# Until now, this motion has typically consisted of a single `Motion` component.
# In this case, it is a `MotionList` containing the contraction (`HeartBeat`) and the disabled,
# zero-amplitude `Rotation` component:

obj.motion

# Now, we will create a bSSFP cine sequence with the following parameters:

RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 32          # image size = N x N
N_phases     = 16          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 6.25e-3     # [s] bSSFP repetition time
frame_interval = RRs[1] / N_phases # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3;     # [s]

# The bSSFP sequence runs continuously throughout each heartbeat. Its ADC is enabled only
# at the 16 evenly spaced cardiac phases, so every k-space line uses the same fixed TR.

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    adc_duration = adc_duration,
);

## Simulation  #hide
raw1 = @suppress simulate(obj, seq, sys); #hide
## Reconstruction #hide
frames1 = @suppress reconstruct_cine(raw1, seq, N_matrix, N_phases); #hide

# The simulation and subsequent reconstruction produces the following cine frames, 
# which look clean and temporally coherent:

fps = 1 / frame_interval #hide
p2 = plot_cine(frames1, fps; Δt=frame_interval); #hide
#md p2 #hide
#jl display(p2);
#nb display(p2);

# ### 2. Arrhythmic Phantom: Variable RR, Constant Sequence
#
# Now, we will introduce arrhythmias into the phantom by varying its RR intervals.
# However, the sequence will still assume a constant RR interval of 1 second.

RRs = [900, 1100, 1000, 1000, 1000, 800] .* 1e-3;

#md # !!! note
#md #     The `RRs` array contains **scaling factors** relative to the original 
#md #     duration of the phantom’s motion cycle.  
#md #     
#md #     In this example, the base duration of the cardiac motion is 1 second, which is defined within the `t` field of its [`TimeCurve`](@ref) structure.   
#md #     Consequently, the elements in `RRs` directly represent the actual RR intervals in seconds (e.g., 0.9 s, 1.1 s, etc.).

# Let's apply the new `RRs` to the phantom:

## Take the time curve from the contraction motion:
t_curve = obj.motion.motions[1].time 
## Generate a new time curve:
t_curve_new = TimeCurve(
    t = t_curve.t, 
    t_unit = t_curve.t_unit,
    periodic = true,
    periods = RRs
)
## Assign the new time curve to both stored motion components:
obj.motion.motions[1].time = t_curve_new
obj.motion.motions[2].time = t_curve_new;

# Let’s visualize how the motion pattern has changed, now with variable-duration RR intervals:

p3 = plot_phantom_map(obj, :T1 ; height=450, time_samples=41); #hide
#md p3 #hide

#jl display(p3);
#nb display(p3);

# Since the sequence still assumes a constant RR interval, it becomes unsynchronized with the phantom.
# This results in artifacts and temporal inconsistencies in the cine images. We will showcase these images in the next section.

## Simulation  #hide
raw2 = @suppress simulate(obj, seq, sys) #hide
## Reconstruction #hide
frames2 = @suppress reconstruct_cine(raw2, seq, N_matrix, N_phases); #hide
#jl display(plot_cine(frames2, fps; Δt=frame_interval));
#nb display(plot_cine(frames2, fps; Δt=frame_interval));

# ### 3. Prospective Triggering: Resynchronized Acquisition
# To correct this, we synchronize the sequence **manually** by providing it the same RR intervals as the phantom:
seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    adc_duration = adc_duration,
);

# This approach manually mimics cardiac triggering. The resulting cine is 
# once again correctly aligned, despite the underlying arrhythmia.
# 
# In the future, this synchronization will be handled automatically 
# through scanner trigger extensions; this tutorial keeps the synchronization explicit.

## Simulation  #hide
raw3 = @suppress simulate(obj, seq, sys) #hide
## Reconstruction #hide
frames3 = @suppress reconstruct_cine(raw3, seq, N_matrix, N_phases); #hide

#jl display(plot_cine(frames3, fps; Δt=frame_interval));
#nb display(plot_cine(frames3, fps; Δt=frame_interval));

#md # Below, we compare the results of the desynchronized 👎 acquisition simulated in the previous section with the resynchronized 🕐 acquisition: 
p4 = plot_cine([frames2 ;; frames3], fps; Δt=frame_interval); #hide
#md p4 #hide
#jl display(p4);
#nb display(p4);

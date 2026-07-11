# # Cardiac Cine MRI with Arrhythmias

using KomaMRI, Plots, Printf, Suppressor; #hide
import KomaMRI.PulseDesigner as PD #hide
sys = Scanner(); #hide

function bSSFP_cine(FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; #hide
    Δf=0, adc_duration=1e-3, N_dummy_cycles=10) #hide
    RR(i) = RRs[(i - 1) % length(RRs) + 1] #hide
    seq = Sequence() #hide
    base_seq = bSSFP(FOV, N_matrix, TR, flip_angle, sys; Δf, adc_duration) #hide
    n = 1 #hide
    for i in 0:(N_matrix - 1) #hide
        line = base_seq[6 * i .+ (1:6)] #hide
        if N_dummy_cycles > 0 && i == 0 #hide
            dummy_dur = N_dummy_cycles * dur(line) #hide
            rr_sum = RR(n) #hide
            while dummy_dur > rr_sum #hide
                n += 1 #hide
                rr_sum += RR(n) #hide
            end #hide
            @addblock seq += Delay(rr_sum - dummy_dur) #hide
        end #hide
        for j in 1:(N_phases + N_dummy_cycles) #hide
            l = copy(line) #hide
            l[1].RF[1].A *= (-1)^j #hide
            l[4].ADC[1].ϕ = iseven(j) ? 0 : π #hide
            if j <= N_dummy_cycles #hide
                l = 0 * l #hide
                l[4].ADC[1].N = 0 #hide
            end #hide
            @addblock seq += l #hide
        end #hide
        phase_dur = (N_phases + N_dummy_cycles) * dur(line) #hide
        n += 1 #hide
        rr_sum = RR(n) #hide
        while phase_dur > rr_sum #hide
            n += 1 #hide
            rr_sum += RR(n) #hide
        end #hide
        @addblock seq += Delay(rr_sum - phase_dur) #hide
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

function plot_cine(frames, fps; Δt=1 / fps, filename="cine_recon.gif", width=400, height=400) #hide
    x = 0:(size(frames[1])[2] - 1) #hide
    y = 1:size(frames[1])[1] #hide
    global_min = minimum(reduce(vcat, frames)) #hide
    global_max = maximum(reduce(vcat, frames)) #hide
    n_frames, n_cines = ndims(frames) == 1 ? (length(frames), 1) : size(frames) #hide
    t = 0 #hide
    anim = @animate for i in 1:n_frames #hide
        t += Δt #hide
        plots = [ #hide
            Plots.plot!( #hide
                Plots.heatmap(x, y, frames[i, j]', color=:greys; aspect_ratio=:equal, colorbar=true, clim=(global_min, global_max), size=(n_cines * width, height)), #hide
                title="t = " * Printf.@sprintf("%.3f", t) * "s", #hide
                xlims=(minimum(x), maximum(x)), #hide
                ylims=(minimum(y), maximum(y)), #hide
            ) for j in 1:n_cines #hide
        ] #hide
        Plots.plot(plots..., layout=(1, n_cines)) #hide
    end #hide
    gif(anim, filename, fps=fps) #hide
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

# This tutorial shows how to simulate cardiac cine MRI using Koma,  
# including cases with variable RR intervals (i.e., arrhythmias). You'll learn how to:

# 1. Simulate a clean cine acquisition with constant RR intervals.  
# 2. Introduce arrhythmias (variable RR intervals) into the cardiac phantom.  
# 3. Observe how this desynchronization degrades image quality.  
# 4. Correct the acquisition by synchronizing the sequence with the phantom’s RR variability (manual triggering).

# ### 1. Constant RR for Phantom and Sequence
#
# We will begin by simulating a cardiac cine on a myocardial phantom with a constant RR interval.
# We'll use the `heart_phantom` function to create a ring-shaped phantom filled with blood, resembling the left ventricle:

obj = heart_phantom(; spins_per_voxel=20);

# By default, this phantom exhibits periodic contraction and rotation with a 1-second period:
p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21) #hide
#jl display(p1);

# As shown in previous tutorials, the phantom's motion is defined by its `motion` field.
# Until now, this motion has typically consisted of a single `Motion` component.
# In this case, however, it consists of two independent motions: a contraction (`HeartBeat`)
# and a `Rotation`. These two are grouped together in a `MotionList` structure:

obj.motion

# Now, we will create a bSSFP cine sequence with the following parameters:

RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 32          # image size = N x N
N_phases     = 16          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 25e-3       # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3;     # [s]

# 

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 16, adc_duration = adc_duration,
);

## Simulation  #hide
raw1 = @suppress simulate(obj, seq, sys); #hide
## Reconstruction #hide
frames1 = @suppress reconstruct_cine(raw1, seq, N_matrix, N_phases); #hide

# The simulation and subsequent reconstruction produces the following cine frames, 
# which look clean and temporally coherent:

fps = 25 #hide
p2 = @suppress plot_cine(frames1, fps; Δt=TR, filename="../public/assets/tut-7-frames1.gif"); #hide
#jl display(p2);
#nb display(p2);

#md # ```@raw html
#md # <center><object data="../assets/tut-7-frames1.gif" style="width:100%; max-width:325px"></object></center>
#md # ```

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
## Assign the new time curve to both the contraction and the rotation:
obj.motion.motions[1].time = t_curve_new
obj.motion.motions[2].time = t_curve_new;

# Let’s visualize how the motion pattern has changed, now with variable-duration RR intervals:

p3 = plot_phantom_map(obj, :T1 ; height=450, time_samples=41) #hide

#jl display(p3);

# Since the sequence still assumes a constant RR interval, it becomes unsynchronized with the phantom.
# This results in artifacts and temporal inconsistencies in the cine images. We will showcase these images in the next section.

## Simulation  #hide
raw2 = @suppress simulate(obj, seq, sys) #hide
## Reconstruction #hide
frames2 = @suppress reconstruct_cine(raw2, seq, N_matrix, N_phases); #hide
@suppress plot_cine(frames2, fps; Δt=TR, filename="../public/assets/tut-7-frames2.gif"); #hide
#jl plot_cine(frames2, fps; Δt=TR, filename="tut-7-frames2.gif");
#nb plot_cine(frames2, fps; Δt=TR, filename="tut-7-frames2.gif");

# ### 3. Prospective Triggering: Resynchronized Acquisition
# To correct this, we synchronize the sequence **manually** by providing it the same RR intervals as the phantom:
seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 16, adc_duration = adc_duration,
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

#jl plot_cine(frames3, fps; Δt=TR, filename="tut-7-frames3.gif");
#nb plot_cine(frames3, fps; Δt=TR, filename="tut-7-frames3.gif");

#md # Below, we compare the results of the desynchronized 👎 acquisition simulated in the previous section with the resynchronized 🕐 acquisition: 
@suppress plot_cine([frames2 ;; frames3], fps; Δt=TR, filename="../public/assets/tut-7-frames_comparison.gif"); #hide
#md # ```@raw html
#md # <center><object data="../assets/tut-7-frames_comparison.gif" style="width:100%"></object></center>
#md # ```

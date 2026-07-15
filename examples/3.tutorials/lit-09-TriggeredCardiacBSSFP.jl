# # Triggered Cardiac bSSFP with Arrhythmia

using KomaMRI #hide
import KomaMRI.PulseDesigner as PD #hide
sys = Scanner(); #hide
sim_params = Dict{String,Any}("gpu"=>false); #hide
#hide
function triggered_bSSFP(FOV, N_matrix, TR, flip_angle, sys; #hide
    Δf=0, adc_duration=1e-3, N_ramp_cycles=16, lines_per_heartbeat=8, #hide
    trigger_delay) #hide
    N_matrix % lines_per_heartbeat == 0 || #hide
        error("`lines_per_heartbeat` must divide `N_matrix`.") #hide
    seq = Sequence() #hide
    base_seq = bSSFP(FOV, N_matrix, TR, flip_angle, sys; Δf, adc_duration) #hide
    blocks_per_line = length(base_seq) ÷ N_matrix #hide
    line(i) = base_seq[blocks_per_line * i .+ (1:blocks_per_line)] #hide
    ramp_line = line(N_matrix ÷ 2) #hide
    ramp_line.GR[1:2, :] .*= 0 #hide
    rf_block = findfirst(i -> is_RF_on(ramp_line[i]), 1:length(ramp_line)) #hide
    adc_block = findfirst(i -> is_ADC_on(ramp_line[i]), 1:length(ramp_line)) #hide
    adc = ramp_line.ADC[adc_block] #hide
    first_adc_offset = N_ramp_cycles * dur(ramp_line) + sum(ramp_line.DUR[1:(adc_block - 1)]) + adc.delay + adc.T / 2 #hide
    post_trigger_delay = trigger_delay - first_adc_offset #hide
    post_trigger_delay >= 0 || error("The flip-angle ramp does not fit before the requested trigger delay.") #hide
    for first_line in 0:lines_per_heartbeat:(N_matrix - 1) #hide
        @addblock seq += PD.make_trigger(:physio1; duration=post_trigger_delay, sys) #hide
        for cycle in 1:N_ramp_cycles #hide
            l = copy(ramp_line) #hide
            l.RF[rf_block].A *= cycle / N_ramp_cycles * (-1)^cycle #hide
            l.ADC[adc_block].N = 0 #hide
            @addblock seq += l #hide
        end #hide
        last_line = first_line + lines_per_heartbeat - 1 #hide
        for (segment_cycle, line_index) in enumerate(first_line:last_line) #hide
            cycle = N_ramp_cycles + segment_cycle #hide
            l = copy(line(line_index)) #hide
            l.RF[rf_block].A *= (-1)^cycle #hide
            l.ADC[adc_block].ϕ = iseven(cycle) ? 0 : π #hide
            @addblock seq += l #hide
        end #hide
    end #hide
    return seq #hide
end #hide
#hide
function bSSFP(FOV, N, TR, flip_angle, sys; #hide
    Δf=0, pulse_duration=2e-3, z0=0.0, slice_thickness=10e-3, TBP=4, adc_duration=1e-3) #hide
    Δk = 1 / FOV #hide
    ro_area = N * Δk #hide
    rf, gz, gz_reph = PD.make_sinc_pulse( #hide
        flip_angle * π / 180; #hide
        duration=pulse_duration, #hide
        slice_thickness, #hide
        freq_offset=Δf + TBP / pulse_duration * z0 / slice_thickness, #hide
        time_bw_product=TBP, #hide
        apodization=0.5, #hide
        sys, #hide
    ) #hide
    readout_time = max(adc_duration, N * sys.ADC_Δt) #hide
    readout_time == adc_duration || @warn "ADC duration is too short. It will be extended to $(readout_time * 1e3) ms." #hide
    gx = PD.make_trapezoid(; flat_area=ro_area, flat_time=readout_time, sys) #hide
    adc = PD.make_adc(N; duration=gx.T, delay=gx.rise, sys) #hide
    excitation_time = max(dur(rf), dur(gz)) #hide
    pre_time = dur(gz_reph) #hide
    gx_pre = PD.make_trapezoid(; area=-area(gx) * γ / 2, duration=pre_time, sys) #hide
    bssfp = Sequence(sys) #hide
    @addblock for i in 0:(N - 1) #hide
        ky = (i - N / 2) * Δk #hide
        gy_pre = PD.make_trapezoid(; area=ky, duration=pre_time, rise_time=gx_pre.rise, fall_time=gx_pre.fall, sys) #hide
        delay_TR = TR - excitation_time - pre_time - dur(gx) - pre_time #hide
        delay_before_readout = round_to_raster( #hide
            TR / 2 - (excitation_time + pre_time + adc.delay + adc.T / 2 - rf.delay - rf.center), #hide
            sys.GR_Δt, #hide
        ) #hide
        delay_after_readout = delay_TR - delay_before_readout #hide
        bssfp += (rf, z=gz) #hide
        bssfp += (x=gx_pre, y=gy_pre, z=gz_reph) #hide
        bssfp += Delay(delay_before_readout) #hide
        bssfp += (adc, x=gx) #hide
        bssfp += Delay(delay_after_readout) #hide
        bssfp += (x=gx_pre, y=-gy_pre, z=gz_reph) #hide
    end #hide
    bssfp.DEF = Dict("Nx"=>N, "Ny"=>N, "Nz"=>1, "Name"=>"bssfp$(N)x$(N)", "FOV"=>[FOV, FOV, 0]) #hide
    return bssfp #hide
end #hide
#hide
function reconstruct_image(raw, seq, N_matrix) #hide
    acq = AcquisitionData(raw) #hide
    _, ktraj = get_kspace(seq) #hide
    acq.traj[1].circular = false #hide
    acq.traj[1].nodes = permutedims(ktraj[:, 1:2]) #hide
    acq.traj[1].nodes ./= maximum(2 * abs.(acq.traj[1].nodes[:])) #hide
    acq.traj[1].numProfiles = N_matrix #hide
    rec_params = Dict{Symbol,Any}( #hide
        :reco=>"direct", :reconSize=>(N_matrix, N_matrix), :densityWeighting=>false, #hide
    ) #hide
    image = reconstruction(acq, rec_params).data #hide
    return abs.(reshape(image, N_matrix, N_matrix, :)[:, :, 1]) #hide
end #hide
#hide
function cardiac_phantom(; RR_intervals) #hide
    base_RR = 1.0 #hide
    systole_end = 0.42 #hide
    rest_start = 0.67 #hide
    rest_end = 0.77 #hide
    rest_position = 1 - (rest_start - systole_end) / (1 - systole_end - (rest_end - rest_start)) #hide
    motion_time = TimeCurve( #hide
        t=base_RR .* [0.0, systole_end, rest_start, rest_end, 1.0], #hide
        t_unit=[0.0, 1.0, rest_position, rest_position, 0.0], #hide
        periodic=true, #hide
        periods=RR_intervals ./ base_RR, #hide
    ) #hide
    obj = heart_phantom(; rotation_angle=5.0, spins_per_voxel=20) #hide
    obj.motion = MotionList(( #hide
        Motion(motion.action, motion_time, motion.spins) for motion in obj.motion.motions #hide
    )...) #hide
    return obj #hide
end #hide
nothing #hide
#hide
# In this tutorial, you will acquire one 32×32 mid-diastolic bSSFP image using
# 8 phase-encoding lines per heartbeat. You will then simulate an arrhythmia while the
# sequence continues to trigger every second and observe the resulting artifact.

# ### 1. Correct triggering with a regular cardiac rhythm
#
# We begin with a myocardial phantom whose contraction and rotation repeat every second.
# It contracts during systole, relaxes during diastole, and pauses around the
# mid-diastolic acquisition window.

obj = cardiac_phantom(; RR_intervals=1.0)
trigger_delay = 0.72;

# The resulting contraction and rotation have a one-second period:
p1 = plot_phantom_map(obj, :T1; height=450, max_spins=2000, time_samples=21) #hide
#jl display(p1);

# Each heartbeat starts with a flip-angle ramp. The first of 8 phase-encoding lines is
# acquired 720 ms after the R peak; the remaining lines follow every TR within the same
# mid-diastolic window.

N_matrix           = 32          # image size = N x N
lines_per_heartbeat = 8
N_ramp_cycles       = 16          # flip-angle ramp length
FOV                 = 0.11        # [m]
TR                  = 2.96e-3     # [s]
flip_angle          = 40          # [º]
adc_duration        = 0.24e-3;    # [s]

seq = triggered_bSSFP(
    FOV, N_matrix, TR, flip_angle, sys;
    N_ramp_cycles, lines_per_heartbeat, trigger_delay, adc_duration,
);

# A heart rate of `1` means 1 Hz, giving one R peak per second:

regular_cardiac_rythm = CardiacSignal(; heart_rate=1);

# Passing the signal to `plot_seq` resolves the trigger waits and adds the ECG trace.
# The complete sequence is plotted, with the initial view showing the first two seconds:

pseq = plot_seq(
    seq;
    physio=regular_cardiac_rythm,
    range=[0, 2000],
    gl=true,
    show_adc=true,
    height=450,
);
pseq #hide

# Simulate the four triggered heartbeat segments:

raw_regular = simulate(
    obj, seq, sys;
    sim_params,
    physio=regular_cardiac_rythm,
    verbose=false,
);

image_regular = reconstruct_image(raw_regular, seq, N_matrix); #hide
pimage_regular = plot_image( #hide
    image_regular; height=400, title="Correctly triggered regular rhythm", #hide
); #hide

# The segments combine into one clean mid-diastolic image:

#md pimage_regular #hide
#jl display(pimage_regular);
#nb display(pimage_regular);

# ### 2. Sequence triggered incorrectly during arrhythmia
#
# Now make the phantom arrhythmic while keeping the trigger signal at 1 Hz:

arrhythmic_intervals = [1.2, 1.2, 0.8, 0.8]
obj_arrhythmic = cardiac_phantom(; RR_intervals=arrhythmic_intervals);

# The varying cycle lengths are visible in the motion plot:

p3 = plot_phantom_map(obj_arrhythmic, :T1; height=450, max_spins=2000, time_samples=41) #hide
#jl display(p3);

# The regular trigger now samples the segments at different cardiac phases:

raw_mismatched = simulate(
    obj_arrhythmic, seq, sys;
    sim_params,
    physio=regular_cardiac_rythm, # triggered incorrectly every second
    verbose=false,
);

image_mismatched = reconstruct_image(raw_mismatched, seq, N_matrix); #hide
image_scale = maximum(image_mismatched) #hide
pimage_mismatched = plot_image( #hide
    image_mismatched; height=400, zmin=0, zmax=image_scale, #hide
    title="Incorrect sequence triggering", #hide
); #hide

# The phase-encoding segments no longer describe the same cardiac state:

#md pimage_mismatched #hide
#jl display(pimage_mismatched);
#nb display(pimage_mismatched);

# The phantom's cardiac rhythm and `CardiacSignal` are independent inputs: changing the
# motion does not alter the trigger schedule.

# # Cardiac Cine MRI with Arrhythmias

using KomaMRI, PlotlyJS, Plots, Printf, Suppressor; #hide
include(joinpath(dirname(pathof(KomaMRI)), "../examples/3.tutorials/utils/RRVariability.jl")); #hide
sys = Scanner(); #hide

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

obj = heart_phantom(); 

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
N_matrix     = 50          # image size = N x N
N_phases     = 25          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 25e-3       # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3;     # [s]

# 

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 40, adc_duration = adc_duration,
);

## Simulation  #hide
raw1 = @suppress simulate(obj, seq, sys); #hide
## Reconstruction #hide
frames1 = @suppress reconstruct_cine(raw1, seq, N_matrix, N_phases); #hide

# The simulation and subsequent reconstruction produces the following cine frames, 
# which look clean and temporally coherent:

fps = 25 #hide
p2 = @suppress plot_cine(frames1, fps; Δt=TR, filename="../assets/tut-7-frames1.gif"); #hide
#jl display(p2);
#nb display(p2);

#md # ```@raw html
#md # <center><object data="../../assets/tut-7-frames1.gif" style="width:100%; max-width:325px"></object></center>
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

#jl plot_cine(frames2, fps; Δt=TR, filename="tut-7-frames2.gif");
#nb plot_cine(frames2, fps; Δt=TR, filename="tut-7-frames2.gif");

# ### 3. Prospective Triggering: Resynchronized Acquisition
# To correct this, we synchronize the sequence **manually** by providing it the same RR intervals as the phantom:
seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 40, adc_duration = adc_duration,
);

# This approach manually mimics cardiac triggering. The resulting cine is 
# once again correctly aligned, despite the underlying arrhythmia.
# 
# In the future, this synchronization will be handled automatically 
# through upcoming support for trigger extensions in the sequence framework.

## Simulation  #hide
raw3 = @suppress simulate(obj, seq, sys) #hide
## Reconstruction #hide
frames3 = @suppress reconstruct_cine(raw3, seq, N_matrix, N_phases); #hide

#jl plot_cine(frames3, fps; Δt=TR, filename="tut-7-frames3.gif");
#nb plot_cine(frames3, fps; Δt=TR, filename="tut-7-frames3.gif");

#md # Below, we compare the results of the desynchronized 👎 acquisition simulated in the previous section with the resynchronized 🕐 acquisition: 
#md @suppress plot_cine([frames2 ;; frames3], fps; Δt=TR, filename="../assets/tut-7-frames_comparison.gif"); #hide
#md # ```@raw html
#md # <center><object data="../../assets/tut-7-frames_comparison.gif" style="width:100%"></object></center>
#md # ```


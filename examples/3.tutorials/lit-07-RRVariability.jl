# # Arrhythmic Cardiac Cine MRI

using KomaMRI, PlotlyJS, Plots, Printf # hide

include(joinpath(dirname(pathof(KomaMRI)), "../examples/3.tutorials/utils/RRVariability.jl")) # hide

sys = Scanner() # hide

# This tutorial shows how to simulate cardiac cine MRI using Koma.  
# including cases with variable RR intervals (i.e., arrhythmias). You'll learn how to:

# 1. Simulate a clean cine acquisition with constant RR intervals.
# 2. Introduce arrhythmias (variable RR) in the cardiac phantom.
# 3. Observe how this desynchronization degrades image quality.
# 4. Fix the acquisition by synchronizing the sequence with the phantom's RR variability (manual triggering).

# ### 1. Constant RR for Phantom and Sequence
# 
# We will begin by simulating a cardiac cine on a myocardial phantom with a constant RR interval.
# Let's call the `heart_phantom` function to create a ring-shaped phantom filled with blood, which resembles the left ventricle:
obj = heart_phantom()

# By default, this phantom exhibits periodic contraction and rotation, with a period of 1 second:
p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21) # hide
#md PlotlyJS.savefig(p1, "../assets/tut-6-phantom.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/tut-6-phantom.html" style="width:90%; height:470px;"></object></center>
#md # ```

# Now, we will create a bSSFP cine sequence with the following parameters:

RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 50          # image size = N x N
N_phases     = 40          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 20e-3       # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3

# 

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 40, adc_duration = adc_duration,
)

## Simulation  # hide
raw1 = simulate(obj, seq, sys) # hide
## Reconstruction # hide
frames1 = reconstruct_cine(raw1, seq, N_matrix, N_phases) # hide

# The simulation and subsequent reconstruction produces the following cine frames, 
# which look clean and temporally coherent:

fps = 25 # hide
p2 = plot_cine(frames1, fps; Δt=TR, filename="../assets/tut-7-frames1.gif"); # hide
#jl #nb display(p2)

#md # ```@raw html
#md # <center><object data="../../assets/tut-7-frames1.gif" style="width:40%"></object></center>
#md # ```

# ### 2. Arrhythmic Phantom: Variable RR, Constant Sequence
# 
# Now, we will introduce arrhythmias in the phantom by modifying its RR intervals. 
# However, we will keep the sequence with a constant RR interval of 1 second:

RRs = [900, 1100, 1000, 1000, 1000, 800] .* 1e-3

#md # !!! note
#md #     The `RRs` array contains **scaling factors** relative to the original 
#md #     duration of the phantom’s motion cycle.  
#md #     
#md #     In this example, the base duration of the cardiac motion is 1 second, which is defined within the `t` field of its [`TimeCurve`](@ref) structure.   
#md #     Consequently, the elements in `RRs` directly represent the actual RR intervals in seconds (e.g., 0.9 s, 1.1 s, etc.).

## Apply the new RRs to the phantom (both contraction and rotation):
obj.motion.motions[1].time.periods = RRs 
obj.motion.motions[2].time.periods = RRs

# Since the sequence still assumes a constant RR, it becomes unsynchronized with the phantom.
# This produces artifacts and temporal inconsistencies in the cine.

## Simulation  # hide
raw2 = simulate(obj, seq, sys) # hide
## Reconstruction # hide
frames2 = reconstruct_cine(raw2, seq, N_matrix, N_phases) # hide

#jl #nb plot_cine(frames2, fps; Δt=TR, filename="tut-7-frames2.gif")

# ### 3. Prospective Triggering: Resynchronized Acquisition
# To correct this, we synchronize the sequence by providing it the same RR variability as the phantom:
seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys; 
    N_dummy_cycles = 40, adc_duration = adc_duration,
)

# This manually mimics cardiac triggering. The resulting cine is 
# once again correctly aligned, despite the underlying arrhythmia.

## Simulation  # hide
raw3 = simulate(obj, seq, sys) # hide
## Reconstruction # hide
frames3 = reconstruct_cine(raw3, seq, N_matrix, N_phases) # hide

#jl #nb plot_cine(frames3, fps; Δt=TR, filename="tut-7-frames3.gif") # hide

#md # Below we compare the results from the desynchronized 👎 and resynchronized 🕐 acquisitions:
#md plot_cine([frames2 ;; frames3], fps; Δt=TR, filename="../assets/tut-7-frames_comparison.gif"); #hide
#md # ```@raw html
#md # <center><object data="../../assets/tut-7-frames_comparison.gif" style="width:80%"></object></center>
#md # ```


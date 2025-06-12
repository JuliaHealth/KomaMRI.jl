using KomaMRI, PlotlyJS, Plots, Printf # hide

include(joinpath(dirname(pathof(KomaMRI)), "../examples/3.tutorials/utils/RRVariability.jl")) # hide

sys = Scanner() # hide

obj = heart_phantom()

p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21) # hide
display(p1)

RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 50          # image size = N x N
N_phases     = 40          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 20e-3       # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    N_dummy_cycles = 40, adc_duration = adc_duration,
)

# Simulation  # hide
raw1 = simulate(obj, seq, sys) # hide
# Reconstruction # hide
frames1 = reconstruct_cine(raw1, seq, N_matrix, N_phases) # hide

fps = 25 # hide
p2 = plot_cine(frames1, fps; Δt=TR, filename="../assets/tut-7-frames1.gif"); # hide

RRs = [900, 1100, 1000, 1000, 1000, 800] .* 1e-3


# Apply the new RRs to the phantom (both contraction and rotation):
obj.motion.motions[1].time.periods = RRs
obj.motion.motions[2].time.periods = RRs

# Simulation  # hide
raw2 = simulate(obj, seq, sys) # hide
# Reconstruction # hide
frames2 = reconstruct_cine(raw2, seq, N_matrix, N_phases) # hide

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    N_dummy_cycles = 40, adc_duration = adc_duration,
)

# Simulation  # hide
raw3 = simulate(obj, seq, sys) # hide
# Reconstruction # hide
frames3 = reconstruct_cine(raw3, seq, N_matrix, N_phases) # hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

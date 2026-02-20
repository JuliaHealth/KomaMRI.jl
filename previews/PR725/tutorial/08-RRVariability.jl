using KomaMRI, PlotlyJS, Plots, Printf, Suppressor; #hide
include(joinpath(dirname(pathof(KomaMRI)), "../examples/3.tutorials/utils/RRVariability.jl")); #hide
sys = Scanner(); #hide

obj = heart_phantom();

p1 = plot_phantom_map(obj, :T1 ; height=450, time_samples=21) #hide
display(p1);

obj.motion

RRs          = [1.0]       # [s] constant RR interval
N_matrix     = 50          # image size = N x N
N_phases     = 25          # Number of cardiac phases
FOV          = 0.11        # [m]
TR           = 25e-3       # [s]
flip_angle   = 10          # [º]
adc_duration = 0.2e-3;     # [s]

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    N_dummy_cycles = 40, adc_duration = adc_duration,
);

# Simulation  #hide
raw1 = @suppress simulate(obj, seq, sys); #hide
# Reconstruction #hide
frames1 = @suppress reconstruct_cine(raw1, seq, N_matrix, N_phases); #hide

fps = 25 #hide
p2 = @suppress plot_cine(frames1, fps; Δt=TR, filename="../assets/tut-7-frames1.gif"); #hide
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
# Assign the new time curve to both the contraction and the rotation:
obj.motion.motions[1].time = t_curve_new
obj.motion.motions[2].time = t_curve_new;

p3 = plot_phantom_map(obj, :T1 ; height=450, time_samples=41) #hide

display(p3);

# Simulation  #hide
raw2 = @suppress simulate(obj, seq, sys) #hide
# Reconstruction #hide
frames2 = @suppress reconstruct_cine(raw2, seq, N_matrix, N_phases); #hide

plot_cine(frames2, fps; Δt=TR, filename="tut-7-frames2.gif");

seq = bSSFP_cine(
    FOV, N_matrix, TR, flip_angle, RRs, N_phases, sys;
    N_dummy_cycles = 40, adc_duration = adc_duration,
);

# Simulation  #hide
raw3 = @suppress simulate(obj, seq, sys) #hide
# Reconstruction #hide
frames3 = @suppress reconstruct_cine(raw3, seq, N_matrix, N_phases); #hide

plot_cine(frames3, fps; Δt=TR, filename="tut-7-frames3.gif");

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

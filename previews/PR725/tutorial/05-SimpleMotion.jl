using KomaMRI #hide
using PlotlyJS, Suppressor #hide
sys = Scanner(); #hide

obj = brain_phantom2D()
obj.Δw .= 0; #hide

obj.motion = translate(2e-2, 0.0, 0.0, TimeRange(t_start=0.0, t_end=200e-3))
p1 = plot_phantom_map(obj, :T2 ; height=450, time_samples=4); #hide
display(p1);
# Read Sequence #hide
seq_file1 = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq") #hide
seq1 = @suppress read_seq(seq_file1); #hide
# Simulate #hide
raw1 = @suppress simulate(obj, seq1, sys) #hide
# Recon #hide
acq1 = AcquisitionData(raw1) #hide
acq1.traj[1].circular = false #hide
Nx, Ny = raw1.params["reconSize"][1:2] #hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) #hide
image1 = reconstruction(acq1, reconParams); #hide

# Plotting the recon #hide
p2 = plot_image(abs.(image1[:, :, 1]); height=400); #hide
display(p2);

sample_times = get_adc_sampling_times(seq1)
displacements = hcat(get_spin_coords(obj.motion, [0.0], [0.0], [0.0], sample_times)...)
p3 = plot( #hide
    sample_times, #hide
    displacements .* 1e2, #hide
    Layout( #hide
        title = "Head displacement in x, y and z", #hide
        xaxis_title = "time (s)", #hide
        yaxis_title = "Displacement (cm)" #hide
    )) #hide
restyle!(p3,1:3, name=["ux(t)", "uy(t)", "uz(t)"]); #hide
display(p3);

_, kspace = get_kspace(seq1)
ΔΦ = 2π*sum(kspace .* displacements, dims=2);

acq1.kdata[1] .*= exp.(im*ΔΦ)
image2 = reconstruction(acq1, reconParams) #hide
p4 = plot_image(abs.(image2[:, :, 1]); height=400) #hide
display(p2)
display(p4);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

using KomaMRI, Suppressor #hide
sys = Scanner(); #hide

obj = brain_phantom2D() # a slice of a brain
p1 = plot_phantom_map(obj, :T2 ; height=400, width=400, view_2d=true)
p2 = plot_phantom_map(obj, :Δw ; height=400, width=400, view_2d=true)
display([p1 p2])

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = @suppress read_seq(seq_file)
p3 = plot_seq(seq; range=[0 40], slider=true, height=300)
display(p3);

raw = @suppress simulate(obj, seq, sys)
p4 = plot_signal(raw; range=[98.4 103.4] , height=300)
display(p4);

# Get the acquisition data
acq = AcquisitionData(raw)
acq.traj[1].circular = false #This is to remove the circular mask

# Setting up the reconstruction parameters
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

# Plotting the recon
slice_abs = abs.(image[:, :, 1])
p5 = plot_image(slice_abs; height=400)
display(p5);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

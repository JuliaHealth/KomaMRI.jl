using KomaMRIPlots
using KomaMRICore

sys = Scanner()

obj = brain_phantom2D()
plot_phantom_map(obj, :T2 ; height=400)

seq_file = "../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq"
seq = read_seq(seq_file)
plot_seq(seq; range=[0 40], slider=true, height=300)
plot_kspace(seq)
plot_M0(seq)

raw = simulate(obj, seq, sys)
plot_signal(raw)

# Get the acquisition data
acq = AcquisitionData(raw)
acq.traj[1].circular = false #This is to remove a circular mask

# Setting up the reconstruction parameters
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

# Plotting the recon
slice_abs = abs.(image[:, :, 1])
p5 = plot_image(slice_abs; height=400)

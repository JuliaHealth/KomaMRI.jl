using KomaMRI # hide
sys = Scanner() # hide

obj = brain_phantom3D()
obj.Δw .= 0 # Removes the off-resonance

obj.motion = SimpleMotion([Rotation(t_start=0.0, t_end=0.5, pitch=15.0, roll=0.0, yaw=45.0)])
p1 = plot_phantom_map(obj, :T2 ; height=600, intermediate_time_samples=4)
display(p1)

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = read_seq(seq_file)
p2 = plot_seq(seq; range=[0 40], slider=true, height=300)
display(p2)

# Simulate
raw1 = simulate(obj, seq, sys)
raw2 = simulate(obj, Delay(0.5) + seq, sys)

# Get the acquisition data
acq1 = AcquisitionData(raw1)
acq2 = AcquisitionData(raw2)
acq1.traj[1].circular = false #This is to remove the circular mask
acq2.traj[1].circular = false

# Setting up the reconstruction parameters
Nx, Ny = raw1.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))

image1 = reconstruction(acq1, reconParams)
image2 = reconstruction(acq2, reconParams)

# Plotting the recon
p3 = plot_image(abs.(image1[:, :, 1]); height=400)
p4 = plot_image(abs.(image2[:, :, 1]); height=400)
display(p3)
display(p4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

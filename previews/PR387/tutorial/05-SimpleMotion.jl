using KomaMRI # hide
sys = Scanner() # hide

obj = brain_phantom2D()
obj.Δw .= 0 # hide

obj.motion = SimpleMotion([
    Rotation(t_start=0.0, t_end=200e-3, yaw=20.0, pitch=0.0, roll=0.0)
    ])
p1 = plot_phantom_map(obj, :T2 ; height=600, motion_samples=4)
display(p1)


# Read Sequence # hide
seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq") # hide
seq = read_seq(seq_file) # hide

# Simulate # hide
raw1 = simulate(obj, seq, sys) # hide

# Recon # hide
acq1 = AcquisitionData(raw1) # hide
acq1.traj[1].circular = false # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

# Plotting the recon
p3 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
display(p3)

# Read Sequence # hide
seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/Spiral/spiral_100x100_FOV230_SPZ_INTER1.seq") # hide
seq = read_seq(seq_file) # hide

# Simulate # hide
raw1 = simulate(obj, seq, sys) # hide

# Recon # hide
acq1 = AcquisitionData(raw1) # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

# Plotting the recon # hide
p4 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
display(p4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

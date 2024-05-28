using KomaMRI # hide
sys = Scanner() # hide

obj = brain_phantom2D()
obj.Δw .= 0 # hide

obj.motion = SimpleMotion([
    Rotation(t_start=0.0, t_end=200e-3, yaw=45.0, pitch=0.0, roll=0.0)
    ])
p1 = plot_phantom_map(obj, :T2 ; height=450, intermediate_time_samples=4) # hide

display(p1)


# Read Sequence # hide
seq_file1 = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq") # hide
seq1 = read_seq(seq_file1) # hide

# Simulate # hide
raw1 = simulate(obj, seq1, sys) # hide

# Recon # hide
acq1 = AcquisitionData(raw1) # hide
acq1.traj[1].circular = false # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

# Plotting the recon # hide
p2 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
display(p2)

obj.motion = SimpleMotion([
    Translation(t_start=0.0, t_end=200e-3, dx=2e-2, dy=0.0, dz=0.0)
])
p3 = plot_phantom_map(obj, :T2 ; height=450, intermediate_time_samples=4) # hide
display(p3)


# Simulate # hide
raw1 = simulate(obj, seq1, sys) # hide

# Recon # hide
acq1 = AcquisitionData(raw1) # hide
acq1.traj[1].circular = false # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

sample_times = get_adc_sampling_times(seq1)

# Since translation is a rigid motion,
# we can obtain the displacements only for one spin,
# as the displacements of the rest will be the same.
displacements = hcat(get_spin_coords(obj.motion, [0.0], [0.0], [0.0], sample_times)...)

p4 = KomaMRIPlots.plot( # hide
    sample_times, # hide
    displacements .* 1e2, # hide
    KomaMRIPlots.Layout( # hide
        title = "Head displacement in x, y and z", # hide
        xaxis_title = "time (s)", # hide
        yaxis_title = "Displacement (cm)" # hide
    )) # hide
KomaMRIPlots.restyle!(p4,1:3, name=["Δx", "Δy", "Δz"]) # hide

display(p4)

# Get k-space
_, kspace = get_kspace(seq1)
# Phase correction: ΔΦcor = 2π(kx*Δx + ky*Δy + kz*Δz)
ΔΦ = 2π*sum(kspace .* displacements, dims=2)
# Apply phase correction
acq2 = copy(acq1)
acq2.kdata[1] .*= exp.(im*ΔΦ)
# Reconstruct
image2 = reconstruction(acq2, reconParams)

p5 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
p6 = plot_image(abs.(image2[:, :, 1]); height=400) # hide

display(p5)
display(p6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

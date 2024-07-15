using KomaMRI, MAT, MRIReco, MRIReco.RegularizedLeastSquares, Interpolations

# Scanner and object
sys = Scanner()
obj = brain_phantom2D()

# Circle parameters
n_coils = 8           # Number of coils (a to h)
radius = 0.1          # Radius of the circle

# Precompute coil center positions
angles = range(0, 2π; length=n_coils + 1)[1:(end - 1)] # Equally spaced angles
coil_centers = [(radius * cos(θ), radius * sin(θ)) for θ in angles]

# Generate Gaussian sensitivity for each coil
coil_sensitivities = []
for (cx, cy) in coil_centers
    coil = exp.(-π * ((obj.x .- cx) .^ 2 / 0.007 .+ (obj.y .- cy) .^ 2 / 0.007))
    push!(coil_sensitivities, coil)
end

# Combine into a single array
coil_sens = hcat(coil_sensitivities...)

# Assign to object
obj.coil_sens = coil_sens
sys.rf_coils = ArbitraryRFCoils(obj.x, obj.y, obj.z, complex.(coil_sens), complex.(ones(size(coil_sens))))

seq_file = joinpath(
    dirname(pathof(KomaMRI)),
    "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq",
)
seq = read_seq(seq_file)
# And simulate:

sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = BlochSimple()
raw = simulate(obj, seq, sys; sim_params)

acq = AcquisitionData(raw) # hide
acq.traj[1].circular = false # hide
Nx, Ny = raw.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco => "direct", :reconSize => (Nx, Ny)) # hide
image = reconstruction(acq, reconParams) # hide
# Determine the number of slices in the 5th dimension
num_slices = size(image, 5)

# Create slices and plots dynamically
slices = [abs.(image[:, :, 1, 1, i, 1]) for i in 1:num_slices]
p = [plot_image(slice; height=400) for slice in slices]

##
# Generate Gaussian sensitivity for each coil
FOV = 230e-3
xq = range(-FOV / 2, FOV / 2; length=Nx)
yq = range(-FOV / 2, FOV / 2; length=Ny)
coil_sensitivities = []
for (cx, cy) in coil_centers
    coil = exp.(-π * ((xq .- cx) .^ 2 / 0.007 .+ (yq' .- cy) .^ 2 / 0.007))
    push!(coil_sensitivities, coil)
end

coil_sens_recon = Float32.(hcat([coil[:] for coil in coil_sensitivities]...))
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:reg] = L2Regularization(1.e-3)       # regularization
params[:iterations] = 40
params[:solver] = CGNR
params[:senseMaps] = reshape(complex.(coil_sens_recon), Nx, Ny, 1, size(coil_sens_recon, 2))

# do reconstruction
Ireco = reconstruction(acq, params)
recon_image = plot_image(abs.(Ireco.data[:, :]))
# Combine all coil sensitivities into a single heatmap
combined_coil_sens = sum(abs.(coil_sens_recon), dims=2)  # Sum along the coil dimension

# Reshape the combined sensitivities back to the 2D grid size
combined_coil_sens_2d = reshape(combined_coil_sens, Nx, Ny)

# Plot the combined coil sensitivities
coil_sens_plot = plot_image(combined_coil_sens_2d)

#[p[1] p[2]; p[3] p[4]; recon_image coil_sens_plot]
[p[1] p[2] p[3] p[4]; p[5] p[6] p[7] p[8]; recon_image recon_image coil_sens_plot coil_sens_plot]

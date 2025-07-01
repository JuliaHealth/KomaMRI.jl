 # # Experimental: Simulating with realistic coils

using KomaMRI # hide
coil_spread = 0.02 # hide
coil_offset = 0.1 # hide
coil1(x,y,z) = exp(-π * ((x + coil_offset) ^ 2 / coil_spread + y ^ 2 / coil_spread))
coil2(x,y,z) = exp(-π * ((x - coil_offset) ^ 2 / coil_spread + y ^ 2 / coil_spread))
coil3(x,y,z) = exp(-π * (x ^ 2 / coil_spread + (y + coil_offset) ^ 2 / coil_spread))
coil4(x,y,z) = exp(-π * (x ^ 2 / coil_spread + (y - coil_offset) ^ 2 / coil_spread))

obj = brain_phantom2D()

#FOR ArbitraryRFCoils
max_x = maximum(abs.(obj.x))
max_y = maximum(abs.(obj.y))
max_z = maximum(abs.(obj.z))
x = range(-max_x, max_x, 11)
y = range(-max_x, max_x, 11)
z = range(-max_x, max_x, 11)
coil_sens1 = [coil1(x,y,z) for x in x, y in y, z in z]
coil_sens2 = [coil2(x,y,z) for x in x, y in y, z in z]
coil_sens3 = [coil3(x,y,z) for x in x, y in y, z in z]
coil_sens4 = [coil4(x,y,z) for x in x, y in y, z in z]
coil_sens = cat(complex.(coil_sens1), complex.(coil_sens2), complex.(coil_sens3), complex.(coil_sens4); dims=4)
sys = Scanner(rf_coils = ArbitraryRFCoils(collect(x), collect(y), collect(z), coil_sens, coil_sens))

#FOR RFCoilsSensDefinedAtPhantomPositions
#coil_sens1 = [coil1(x,y,z) for (x,y,z) in zip(obj.x, obj.y, obj.z)]
#coil_sens2 = [coil2(x,y,z) for (x,y,z) in zip(obj.x, obj.y, obj.z)]
#coil_sens3 = [coil3(x,y,z) for (x,y,z) in zip(obj.x, obj.y, obj.z)]
#coil_sens4 = [coil4(x,y,z) for (x,y,z) in zip(obj.x, obj.y, obj.z)]
#coil_sens = cat(complex.(coil_sens1), complex.(coil_sens2), complex.(coil_sens3), complex.(coil_sens4); dims=2)
#sys = Scanner(rf_coils = RFCoilsSensDefinedAtPhantomPositions(coil_sens))


seq_file = joinpath(
    dirname(pathof(KomaMRI)),
    "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq",
)
seq = read_seq(seq_file)

# And simulate:
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = Bloch()
raw = simulate(obj, seq, sys; sim_params)

acq = AcquisitionData(raw) # hide
acq.traj[1].circular = false # hide
Nx, Ny = raw.params["reconSize"][1:2] # hide
# reco parameters
params = Dict{Symbol,Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:iterations] = 40

# Coil sensitivities interpolated to reconstruction points
FOV = 230e-3
xq = range(-FOV / 2, FOV / 2, Nx)
yq = range(-FOV / 2, FOV / 2, Ny)
coil_sens1 = [coil1(x,y,0) for x in xq, y in yq]
coil_sens2 = [coil2(x,y,0) for x in xq, y in yq]
coil_sens3 = [coil3(x,y,0) for x in xq, y in yq]
coil_sens4 = [coil4(x,y,0) for x in xq, y in yq]
coil_sens_recon = Float32.([coil_sens1[:] coil_sens2[:] coil_sens3[:] coil_sens4[:]])
params[:senseMaps] = reshape(complex.(coil_sens_recon), Nx, Ny, 1, size(coil_sens_recon, 2))

# do reconstruction
Ireco = reconstruction(acq, params)
recon_image = plot_image(abs.(Ireco.data[:, :]))
reconParams = Dict{Symbol,Any}(:reco => "direct", :reconSize => (Nx, Ny)) # hide
image = reconstruction(acq, reconParams) # hide
slice_abs1 = abs.(image[:, :, 1, 1, 1, 1]) # hide
slice_abs2 = abs.(image[:, :, 1, 1, 2, 1]) # hide
slice_abs3 = abs.(image[:, :, 1, 1, 3, 1]) # hide
slice_abs4 = abs.(image[:, :, 1, 1, 4, 1]) # hide
p3 = plot_image(slice_abs1; height=400)
p4 = plot_image(slice_abs2; height=400) # hide
p5 = plot_image(slice_abs3; height=400)
p6 = plot_image(slice_abs4; height=400) # hide
[recon_image; p3 p4 p5 p6]
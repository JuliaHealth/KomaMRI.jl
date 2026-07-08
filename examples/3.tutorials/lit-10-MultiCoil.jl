# # Multi-Coil Receive Reconstruction

using KomaMRI, Suppressor #hide

# KomaMRI can model a receive array directly through the [`Scanner`](@ref) object.
# In this tutorial we simulate a 4-coil EPI acquisition, reconstruct the
# individual coil images, and inspect the magnitude and phase of the receive
# sensitivities.

sys = Scanner(rf_rx=BirdcageCoilSens(ncoils=4, radius=0.20, L=0.30))
obj = brain_phantom2D()

seq_file = joinpath(
    dirname(pathof(KomaMRI)),
    "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq",
)
seq = @suppress read_seq(seq_file)

# With a multi-coil receive model attached to `sys`, [`BlochSimple`](@ref)
# acquires one signal per coil at every ADC sample.
sim_params = KomaMRICore.default_sim_params()
sim_params["gpu"] = false
sim_params["sim_method"] = KomaMRICore.BlochSimple()
raw = @suppress simulate(obj, seq, sys; sim_params)

# We reconstruct the coil images directly from the simulated raw data.
acq = AcquisitionData(raw)
acq.traj[1].circular = false
Nx, Ny = raw.params["reconSize"][1:2]
recon_params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, recon_params)

coil_images = image[:, :, 1, 1, :, 1]
coil_abs_images = [abs.(coil_images[:, :, coil]) for coil in 1:size(coil_images, 3)]
rss_image = dropdims(sqrt.(sum(abs2, coil_images; dims=3)); dims=3)
recon_max = maximum(maximum, coil_abs_images) #hide

p1 = plot_image(rss_image; height=400, title="4-coil RSS reconstruction")
#jl display(p1);

# The individual coil reconstructions can be displayed side by side.
p2a = plot_image(coil_abs_images[1]; height=320, title="Coil 1 reconstruction", zmin=0, zmax=recon_max) #hide
p2b = plot_image(coil_abs_images[2]; height=320, title="Coil 2 reconstruction", zmin=0, zmax=recon_max) #hide
p2c = plot_image(coil_abs_images[3]; height=320, title="Coil 3 reconstruction", zmin=0, zmax=recon_max) #hide
p2d = plot_image(coil_abs_images[4]; height=320, title="Coil 4 reconstruction", zmin=0, zmax=recon_max) #hide

# The receive sensitivities are complex-valued. Their magnitude gives the local
# receive weighting and their phase is obtained with `angle.(...)`.
FOVx, FOVy = raw.params["reconFOV"][1:2] .* 1f-3
xs = range(-FOVx / 2, FOVx / 2; length=Nx)
ys = range(-FOVy / 2, FOVy / 2; length=Ny)
xgrid = vec([x for x in xs, _ in ys])
ygrid = vec([y for _ in xs, y in ys])
zgrid = zeros(Float32, length(xgrid))

coil_sens = sys.rf_sens(xgrid, ygrid, zgrid)
coil_sens_abs = [reshape(abs.(coil_sens[coil, :]), Nx, Ny) for coil in 1:size(coil_sens, 1)]
coil_sens_phase = [reshape(angle.(coil_sens[coil, :]), Nx, Ny) for coil in 1:size(coil_sens, 1)]
sens_max = maximum(maximum, coil_sens_abs) #hide

function circularize(img) #hide
    Nx, Ny = size(img)
    cx = (Nx + 1) / 2
    cy = (Ny + 1) / 2
    r = min(Nx, Ny) / 2
    masked = fill(Float32(NaN), Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        if ((i - cx) / r)^2 + ((j - cy) / r)^2 <= 1
            masked[i, j] = img[i, j]
        end
    end
    masked
end #hide

coil_sens_phase_circ = circularize.(coil_sens_phase) #hide
function row_plot(plots, Nx, Ny) #hide
    p = hcat(plots...)
    foreach(t -> t.fields[:showscale] = false, p.plot.data)
    for (i, xref) in enumerate(("x", "x2", "x3", "x4"))
        xaxis = Symbol("xaxis", i)
        yaxis = Symbol("yaxis", i)
        p.plot.layout.fields[yaxis][:scaleanchor] = xref
        p.plot.layout.fields[yaxis][:constrain] = "domain"
        p.plot.layout.fields[xaxis][:range] = [-0.5, Nx - 0.5]
        p.plot.layout.fields[yaxis][:range] = [-0.5, Ny - 0.5]
    end
    p
end #hide

# The sensitivity magnitudes are shown as square maps.
p3a = plot_image(coil_sens_abs[1]; height=320, title="Coil 1 sensitivity magnitude", zmin=0, zmax=sens_max) #hide
p3b = plot_image(coil_sens_abs[2]; height=320, title="Coil 2 sensitivity magnitude", zmin=0, zmax=sens_max) #hide
p3c = plot_image(coil_sens_abs[3]; height=320, title="Coil 3 sensitivity magnitude", zmin=0, zmax=sens_max) #hide
p3d = plot_image(coil_sens_abs[4]; height=320, title="Coil 4 sensitivity magnitude", zmin=0, zmax=sens_max) #hide

# The phase maps keep the circular mask.
p4a = plot_image(coil_sens_phase_circ[1]; height=320, title="Coil 1 sensitivity phase", zmin=-π, zmax=π, colorscale="Jet") #hide
p4b = plot_image(coil_sens_phase_circ[2]; height=320, title="Coil 2 sensitivity phase", zmin=-π, zmax=π, colorscale="Jet") #hide
p4c = plot_image(coil_sens_phase_circ[3]; height=320, title="Coil 3 sensitivity phase", zmin=-π, zmax=π, colorscale="Jet") #hide
p4d = plot_image(coil_sens_phase_circ[4]; height=320, title="Coil 4 sensitivity phase", zmin=-π, zmax=π, colorscale="Jet") #hide

# Putting the coil images together with the sensitivity magnitude and phase
# makes the receive profile of each channel easy to compare.
p2 = row_plot([p2a, p2b, p2c, p2d], Nx, Ny) #hide
p3 = row_plot([p3a, p3b, p3c, p3d], Nx, Ny) #hide
p4 = row_plot([p4a, p4b, p4c, p4d], Nx, Ny) #hide
#md p2 #hide
#jl display(p2);
#md p3 #hide
#jl display(p3);
#md p4 #hide
#jl display(p4);

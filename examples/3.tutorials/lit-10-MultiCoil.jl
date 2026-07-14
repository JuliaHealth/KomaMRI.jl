# # Multi-Coil Receive Reconstruction

using KomaMRI, PlotlyJS, Suppressor #hide

# KomaMRI can attach different receive sensitivity models to the scanner. Here
# we simulate four receive channels with `Bloch()`, first using an idealized
# birdcage model and then spatial sensitivity maps supplied by the user.

obj = brain_phantom2D()

seq_file = joinpath(
    dirname(pathof(KomaMRI)),
    "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq",
)
seq = @suppress read_seq(seq_file)

birdcage = BirdcageCoilSens(ncoils=4, radius=0.20, L=0.30)

coil_spread = 0.08f0 #hide
coil_offset = 0.20f0 #hide
coil1(x, y, z) = exp(-Float32(π) * (((x - coil_offset)^2 + y^2 + z^2) / coil_spread)) *
    cis(angle(complex(-y, -(x - coil_offset)))) #hide
coil2(x, y, z) = exp(-Float32(π) * ((x^2 + (y - coil_offset)^2 + z^2) / coil_spread)) *
    cis(angle(complex(-(y - coil_offset), -x))) #hide
coil3(x, y, z) = exp(-Float32(π) * (((x + coil_offset)^2 + y^2 + z^2) / coil_spread)) *
    cis(angle(complex(-y, -(x + coil_offset)))) #hide
coil4(x, y, z) = exp(-Float32(π) * ((x^2 + (y + coil_offset)^2 + z^2) / coil_spread)) *
    cis(angle(complex(-(y + coil_offset), -x))) #hide

coords = LinRange(-0.23f0, 0.23f0, 17) #hide
zcoords = LinRange(-0.01f0, 0.01f0, 3) #hide
coil_sens = cat( #hide
    [coil1(x, y, z) for x in coords, y in coords, z in zcoords],
    [coil2(x, y, z) for x in coords, y in coords, z in zcoords],
    [coil3(x, y, z) for x in coords, y in coords, z in zcoords],
    [coil4(x, y, z) for x in coords, y in coords, z in zcoords];
    dims=4,
) #hide
arbitrary = ArbitraryCoilSens(coords, coords, zcoords, coil_sens)

function sensitivity_maps(receiver, xgrid, ygrid, zgrid, Nx, Ny) #hide
    sens = get_sens(receiver, xgrid, ygrid, zgrid)
    sens_abs = [reshape(abs.(sens[:, coil]), Nx, Ny) for coil in axes(sens, 2)]
    sens_phase = [reshape(angle.(sens[:, coil]), Nx, Ny) for coil in axes(sens, 2)]
    return sens_abs, sens_phase
end #hide

function case_plot(obj, seq, receiver) #hide
    sys = Scanner(receiver=receiver)
    sim_params = KomaMRICore.default_sim_params()
    sim_params["gpu"] = false
    sim_params["sim_method"] = Bloch()
    raw = @suppress simulate(obj, seq, sys; sim_params, verbose=false)

    acq = AcquisitionData(raw)
    acq.traj[1].circular = false
    Nx, Ny = raw.params["reconSize"][1:2]
    recon_params = Dict{Symbol, Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
    image = reconstruction(acq, recon_params)

    coil_images = image[:, :, 1, 1, :, 1]
    coil_abs_images = [abs.(coil_images[:, :, coil]) for coil in axes(coil_images, 3)]
    recon_min = minimum(minimum, coil_abs_images)
    recon_max = maximum(maximum, coil_abs_images)

    FOVx, FOVy = raw.params["reconFOV"][1:2] .* 1f-3
    xs = range(-FOVx / 2, FOVx / 2; length=Nx)
    ys = range(-FOVy / 2, FOVy / 2; length=Ny)
    xgrid = vec([x for x in xs, _ in ys])
    ygrid = vec([y for _ in xs, y in ys])
    zgrid = zeros(Float32, length(xgrid))

    sens_abs, sens_phase = sensitivity_maps(receiver, xgrid, ygrid, zgrid, Nx, Ny)
    sens_min = minimum(minimum, sens_abs)
    sens_max = maximum(maximum, sens_abs)
    ncoils = length(coil_abs_images)

    plots = (
        [plot_image(image; zmin=recon_min, zmax=recon_max) for image in coil_abs_images],
        [plot_image(image; zmin=sens_min, zmax=sens_max) for image in sens_abs],
        [plot_image(image; zmin=-π, zmax=π, colorscale="Jet") for image in sens_phase],
    )
    p = make_subplots(
        rows=3,
        cols=ncoils,
        subplot_titles=reshape(
            [
                "Coil $(coil) $(kind)"
                for kind in ("reco", "magnitude", "phase") for coil in 1:ncoils
            ],
            1,
            :,
        ),
        horizontal_spacing=0.03,
        vertical_spacing=0.08,
    )
    for (row, panels) in pairs(plots), (col, panel) in pairs(panels)
        for trace in panel.plot.data
            trace.fields[:showscale] = false
            add_trace!(p, trace; row, col)
        end
        index = (row - 1) * ncoils + col
        suffix = index == 1 ? "" : string(index)
        xaxis = Symbol("xaxis", suffix)
        yaxis = Symbol("yaxis", suffix)
        xref = index == 1 ? "x" : "x$(index)"
        Nx, Ny = size(plots[row][col].plot.data[1].fields[:z])
        p.plot.layout[xaxis][:range] = [-0.5, Nx - 0.5]
        p.plot.layout[yaxis][:range] = [-0.5, Ny - 0.5]
        p.plot.layout[yaxis][:scaleanchor] = xref
        p.plot.layout[yaxis][:constrain] = "domain"
    end
    relayout!(p; width="100%", height=900, showlegend=false)
    return p
end #hide

# Each case below shows the same 4 channels as three rows: reconstructed coil
# image, sensitivity magnitude, and sensitivity phase.

# ## Birdcage coil sensitivities

# `BirdcageCoilSens` models coils distributed uniformly around the object. Each
# channel weights the transverse magnetization by a different complex spatial
# sensitivity before the received signals are reconstructed.

p1 = case_plot(obj, seq, birdcage) #hide
#md p1 #hide

# ## Arbitrary coil sensitivities

# `ArbitraryCoilSens` uses sampled complex sensitivity maps instead. Their
# magnitude controls spatial receive strength, while their phase changes the
# phase contributed by each location to a channel.

p2 = case_plot(obj, seq, arbitrary) #hide
#md p2 #hide

# # Multi-Coil Receive Reconstruction

using KomaMRI, Suppressor #hide

# KomaMRI can attach different receive models to the [`Scanner`](@ref). Here we
# compare 4-channel birdcage and arbitrary sensitivities with both
# [`BlochSimple`](@ref) and the CPU implementation of [`Bloch()`](@ref).

obj = brain_phantom2D()

seq_file = joinpath(
    dirname(pathof(KomaMRI)),
    "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq",
)
seq = @suppress read_seq(seq_file)

birdcage = BirdcageCoilSens(ncoils=4, radius=0.20, L=0.30)

coil_spread = 0.02f0 #hide
coil_offset = 0.1f0 #hide
coil1(x, y, z) = exp(-Float32(π) * (((x + coil_offset)^2 + y^2) / coil_spread)) * cis(Float32(π) * x / coil_offset) #hide
coil2(x, y, z) = exp(-Float32(π) * (((x - coil_offset)^2 + y^2) / coil_spread)) * cis(-Float32(π) * x / coil_offset) #hide
coil3(x, y, z) = exp(-Float32(π) * ((x^2 + (y + coil_offset)^2) / coil_spread)) * cis(Float32(π) * y / coil_offset) #hide
coil4(x, y, z) = exp(-Float32(π) * ((x^2 + (y - coil_offset)^2) / coil_spread)) * cis(-Float32(π) * y / coil_offset) #hide

coords = LinRange(-0.12f0, 0.12f0, 17) #hide
zcoords = LinRange(-0.01f0, 0.01f0, 3) #hide
coil_sens = cat( #hide
    [coil1(x, y, z) for x in coords, y in coords, z in zcoords],
    [coil2(x, y, z) for x in coords, y in coords, z in zcoords],
    [coil3(x, y, z) for x in coords, y in coords, z in zcoords],
    [coil4(x, y, z) for x in coords, y in coords, z in zcoords];
    dims=4,
) #hide
arbitrary = KomaMRIBase.ArbitraryRFRxCoils(coords, coords, zcoords, coil_sens, 0f0 + 0f0im)

function row_plot(images, titles, zmin, zmax; colorscale=nothing) #hide
    plots = [
        isnothing(colorscale) ?
            plot_image(images[i]; height=320, title=titles[i], zmin=zmin, zmax=zmax) :
            plot_image(images[i]; height=320, title=titles[i], zmin=zmin, zmax=zmax, colorscale=colorscale)
        for i in eachindex(images)
    ]
    p = hcat(plots...)
    foreach(trace -> trace.fields[:showscale] = false, p.plot.data)
    Nx, Ny = size(images[1])
    for (i, xref) in enumerate(("x", "x2", "x3", "x4"))
        xaxis = Symbol("xaxis", i)
        yaxis = Symbol("yaxis", i)
        p.plot.layout.fields[yaxis][:scaleanchor] = xref
        p.plot.layout.fields[yaxis][:constrain] = "domain"
        p.plot.layout.fields[xaxis][:range] = [-0.5, Nx - 0.5]
        p.plot.layout.fields[yaxis][:range] = [-0.5, Ny - 0.5]
    end
    return p
end #hide

function sensitivity_maps(receiver::BirdcageCoilSens, xgrid, ygrid, zgrid, Nx, Ny) #hide
    sens = get_sens(receiver, xgrid, ygrid, zgrid)
    sens_abs = [reshape(abs.(sens[:, coil]), Nx, Ny) for coil in axes(sens, 2)]
    sens_phase = [reshape(angle.(sens[:, coil]), Nx, Ny) for coil in axes(sens, 2)]
    return sens_abs, sens_phase
end #hide

function sensitivity_maps(receiver::KomaMRIBase.ArbitraryRFRxCoils, xgrid, ygrid, zgrid, Nx, Ny) #hide
    sens = [
        let
            itp = KomaMRIBase.extrapolate(
                KomaMRIBase.GriddedInterpolation(
                    (receiver.x, receiver.y, receiver.z),
                    receiver.coil_sens[:, :, :, coil],
                    KomaMRIBase.Gridded(KomaMRIBase.Linear()),
                ),
                zero(eltype(receiver.coil_sens)),
            )
            reshape(itp.(xgrid, ygrid, zgrid), Nx, Ny)
        end
        for coil in 1:get_n_coils(receiver)
    ]
    return [abs.(coil) for coil in sens], [angle.(coil) for coil in sens]
end #hide

function case_plots(obj, seq, receiver, sim_method) #hide
    sys = Scanner(receiver=receiver)
    sim_params = KomaMRICore.default_sim_params()
    sim_params["gpu"] = false
    sim_params["sim_method"] = sim_method
    raw = @suppress simulate(obj, seq, sys; sim_params, verbose=false)

    acq = AcquisitionData(raw)
    acq.traj[1].circular = false
    Nx, Ny = raw.params["reconSize"][1:2]
    recon_params = Dict{Symbol, Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
    image = reconstruction(acq, recon_params)

    coil_images = image[:, :, 1, 1, :, 1]
    coil_abs_images = [abs.(coil_images[:, :, coil]) for coil in axes(coil_images, 3)]
    recon_max = maximum(maximum, coil_abs_images)

    FOVx, FOVy = raw.params["reconFOV"][1:2] .* 1f-3
    xs = range(-FOVx / 2, FOVx / 2; length=Nx)
    ys = range(-FOVy / 2, FOVy / 2; length=Ny)
    xgrid = vec([x for x in xs, _ in ys])
    ygrid = vec([y for _ in xs, y in ys])
    zgrid = zeros(Float32, length(xgrid))

    sens_abs, sens_phase = sensitivity_maps(receiver, xgrid, ygrid, zgrid, Nx, Ny)
    sens_max = maximum(maximum, sens_abs)
    ncoils = length(coil_abs_images)

    p_phase = row_plot(
        sens_phase,
        ["Coil $(coil) phase" for coil in 1:ncoils],
        -π,
        π;
        colorscale="Jet",
    )
    p_magnitude = row_plot(
        sens_abs,
        ["Coil $(coil) magnitude" for coil in 1:ncoils],
        0,
        sens_max,
    )
    p_image = row_plot(
        coil_abs_images,
        ["Coil $(coil) reconstruction" for coil in 1:ncoils],
        0,
        recon_max,
    )

    return p_phase, p_magnitude, p_image
end #hide

# Each case below shows the same 4 channels as three rows: phase, magnitude,
# and reconstructed coil image.

# ## BlochSimple with BirdcageCoilSens

p1, p2, p3 = case_plots(obj, seq, birdcage, KomaMRICore.BlochSimple()) #hide
#md p1 #hide
#md p2 #hide
#md p3 #hide

# ## BlochSimple with ArbitraryRFRxCoils

p4, p5, p6 = case_plots(obj, seq, arbitrary, KomaMRICore.BlochSimple()) #hide
#md p4 #hide
#md p5 #hide
#md p6 #hide

# ## Bloch CPU with BirdcageCoilSens

p7, p8, p9 = case_plots(obj, seq, birdcage, KomaMRICore.Bloch()) #hide
#md p7 #hide
#md p8 #hide
#md p9 #hide

# ## Bloch CPU with ArbitraryRFRxCoils

p10, p11, p12 = case_plots(obj, seq, arbitrary, KomaMRICore.Bloch()) #hide
#md p10 #hide
#md p11 #hide
#md p12 #hide

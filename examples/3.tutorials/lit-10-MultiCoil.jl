# # Multi-Coil Receive Reconstruction

using KomaMRI, PlotlyJS, Suppressor #hide

# KomaMRI can attach different receive models to the scanner. Here we compare
# 4-channel birdcage and arbitrary sensitivities with both `BlochSimple()` and
# the CPU implementation `Bloch()`.

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

# ## Forward-Model Benchmark

# The next four plots show the same CPU forward-model workload with:
# `brain_phantom3D()[1:10000]`, `PulseDesigner.EPI_example()`, and
# `return_type="mat"`. The hover labels expose the measured median times.

const BENCHMARK_COILS = [1, 2, 4, 8, 16, 32, 64] #hide

function benchmark_plot(title, receiver_label, current_uniform, current_series, upstream_uniform) #hide
    traces = [
        scatter(
            x=[0],
            y=[current_uniform],
            mode="markers",
            name="current uniform",
            marker=attr(color="#4c78a8", size=11, symbol="circle"),
            hovertemplate="current uniform<br>time=%{y:.3f} ms<extra></extra>",
        ),
        scatter(
            x=BENCHMARK_COILS,
            y=current_series,
            mode="lines+markers",
            name="current $(receiver_label)",
            line=attr(color="#1f77b4", width=3),
            marker=attr(size=9),
            hovertemplate="$(receiver_label) %{x} coils<br>time=%{y:.3f} ms<extra></extra>",
        ),
        scatter(
            x=[0],
            y=[upstream_uniform],
            mode="markers",
            name="upstream uniform",
            marker=attr(color="#111111", size=12, symbol="diamond"),
            hovertemplate="upstream uniform<br>time=%{y:.3f} ms<extra></extra>",
        ),
    ]

    return plot(
        traces,
        Layout(
            title=title,
            template="plotly_white",
            width=950,
            height=500,
            margin=attr(l=50, r=50, t=60, b=50),
            legend=attr(orientation="h", y=1.12, x=0.0),
            xaxis=attr(
                title="coil number",
                tickmode="array",
                tickvals=[0, BENCHMARK_COILS...],
                ticktext=["uniform", string.(BENCHMARK_COILS)...],
                range=[-0.6, 66.5],
            ),
            yaxis=attr(title="median time (ms)"),
        ),
        config=PlotConfig(scrollZoom=true),
    )
end #hide

blochsimple_birdcage_bench = benchmark_plot( #hide
    "BlochSimple with birdcage coils",
    "birdcage",
    208.587,
    [247.479, 267.157, 272.183, 362.821, 552.686, 910.303, 1661.54],
    264.729958,
) #hide
#md blochsimple_birdcage_bench #hide

blochsimple_arbitrary_bench = benchmark_plot( #hide
    "BlochSimple with arbitrary coils",
    "arbitrary",
    233.129,
    [271.328, 288.998, 467.817, 570.708, 945.753, 1666.16, 3148.76],
    264.729958,
) #hide
#md blochsimple_arbitrary_bench #hide

blochcpu_birdcage_bench = benchmark_plot( #hide
    "BlochCPU with birdcage coils",
    "birdcage",
    99.7675,
    [322.536, 517.122, 857.677, 2139.26, 6360.35, 12275.4, 23411.8],
    137.578292,
) #hide
#md blochcpu_birdcage_bench #hide

blochcpu_arbitrary_bench = benchmark_plot( #hide
    "BlochCPU with arbitrary coils",
    "arbitrary",
    102.2,
    [419.54, 726.682, 1377.43, 2827.08, 5347.42, 11087.9, 22788.7],
    137.578292,
) #hide
#md blochcpu_arbitrary_bench #hide

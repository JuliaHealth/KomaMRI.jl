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
arbitrary = ArbitraryCoilSens(coords, coords, zcoords, coil_sens)

function sensitivity_maps(receiver::BirdcageCoilSens, xgrid, ygrid, zgrid, Nx, Ny) #hide
    sens = get_sens(receiver, xgrid, ygrid, zgrid)
    sens_abs = [reshape(abs.(sens[:, coil]), Nx, Ny) for coil in axes(sens, 2)]
    sens_phase = [reshape(angle.(sens[:, coil]), Nx, Ny) for coil in axes(sens, 2)]
    return sens_abs, sens_phase
end #hide

function sensitivity_maps(receiver::ArbitraryCoilSens, xgrid, ygrid, zgrid, Nx, Ny) #hide
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

function case_plot(obj, seq, receiver, sim_method) #hide
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
                for kind in ("reconstruction", "magnitude", "phase") for coil in 1:ncoils
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

# ## BlochSimple with BirdcageCoilSens

p1 = case_plot(obj, seq, birdcage, KomaMRICore.BlochSimple()) #hide
#md p1 #hide

# ## BlochSimple with ArbitraryCoilSens

p2 = case_plot(obj, seq, arbitrary, KomaMRICore.BlochSimple()) #hide
#md p2 #hide

# ## Bloch CPU with BirdcageCoilSens

p3 = case_plot(obj, seq, birdcage, KomaMRICore.Bloch()) #hide
#md p3 #hide

# ## Bloch CPU with ArbitraryCoilSens

p4 = case_plot(obj, seq, arbitrary, KomaMRICore.Bloch()) #hide
#md p4 #hide

# ## Forward-Model Benchmark

# These benchmarks use `brain_phantom3D(ss=3)[1:100_000]`,
# `PulseDesigner.EPI_example()`, and Float32 precision. Inputs were constructed
# before timing, each case was warmed up, and the plots report the median of
# three `BenchmarkTools` samples with one evaluation per sample.
#
# **CPU cases:** 10 Julia threads.
#
# **GPU cases:** `Nthreads = 1` on the CPU side, with Metal synchronized after
# every simulation.
#
# **Metal kernels:** 256 GPU work-items per workgroup for precession and
# excitation.
#
# **100,000 spins:** 391 workgroups and 100,096 padded GPU work-items per
# spin-parallel kernel launch. Metal controls the number executed concurrently.
#
# The master baseline is uniform-coil simulation at `ed1d40f`; the branch
# curves include uniform, 1, 4, and 8 receive channels.

const BENCHMARK_COILS = (1, 4, 8) #hide
const BENCHMARK_SPINS = 100_000 #hide
const BENCHMARK_GROUPSIZE = 256 #hide
const BENCHMARK_RADIUS = 0.20 #hide
const BENCHMARK_LENGTH = 0.30 #hide
const BENCHMARK_GRID_RADIUS = 0.23f0 #hide
const BENCHMARK_GRID_POINTS = 17 #hide
const BENCHMARK_COIL_SPREAD = 0.08f0 #hide

const BENCHMARK_METHODS = ( #hide
    (method=BlochSimple(), gpu=Val(false)), #hide
    (method=Bloch(), gpu=Val(false)), #hide
    (method=Bloch(), gpu=Val(true)), #hide
    (method=BlochMagnus1(), gpu=Val(false)), #hide
    (method=BlochMagnus1(), gpu=Val(true)), #hide
    (method=BlochMagnus2(), gpu=Val(false)), #hide
    (method=BlochMagnus2(), gpu=Val(true)), #hide
    (method=BlochMagnus4(), gpu=Val(false)), #hide
    (method=BlochMagnus4(), gpu=Val(true)), #hide
    (method=BlochMagnus6(), gpu=Val(false)), #hide
    (method=BlochMagnus6(), gpu=Val(true)), #hide
) #hide

benchmark_gpu(::Val{false}) = false #hide
benchmark_gpu(::Val{true}) = true #hide
benchmark_synchronize(::Val{false}) = nothing #hide
benchmark_synchronize(::Val{true}) = Metal.synchronize() #hide

benchmark_receiver(::Val{:birdcage}, ncoils) = BirdcageCoilSens( #hide
    ncoils=ncoils, #hide
    radius=BENCHMARK_RADIUS, #hide
    L=BENCHMARK_LENGTH, #hide
) #hide

function benchmark_sensitivity(x, y, z, coil, ncoils) #hide
    θ = 2f0 * Float32(π) * Float32(coil - 1) / Float32(ncoils) #hide
    xc = Float32(BENCHMARK_RADIUS) * cos(θ) #hide
    yc = Float32(BENCHMARK_RADIUS) * sin(θ) #hide
    magnitude = exp( #hide
        -Float32(π) * (((x - xc)^2 + (y - yc)^2 + z^2) / BENCHMARK_COIL_SPREAD), #hide
    ) #hide
    return magnitude * cis(angle(complex(-(y - yc), -(x - xc)))) #hide
end #hide

function benchmark_receiver(::Val{:arbitrary}, ncoils) #hide
    coords = collect( #hide
        LinRange(-BENCHMARK_GRID_RADIUS, BENCHMARK_GRID_RADIUS, BENCHMARK_GRID_POINTS), #hide
    ) #hide
    zcoords = collect(LinRange(-0.01f0, 0.01f0, 3)) #hide
    sensitivities = ComplexF32[ #hide
        benchmark_sensitivity(x, y, z, coil, ncoils) #hide
        for x in coords, y in coords, z in zcoords, coil in 1:ncoils #hide
    ] #hide
    return ArbitraryCoilSens(coords, coords, zcoords, sensitivities) #hide
end #hide

function benchmark_parameters(sim_method, gpu) #hide
    params = KomaMRICore.default_sim_params() #hide
    params["gpu"] = benchmark_gpu(gpu) #hide
    params["Nthreads"] = Threads.nthreads() #hide
    params["precision"] = "f32" #hide
    params["return_type"] = "mat" #hide
    params["sim_method"] = sim_method #hide
    params["gpu_groupsize_precession"] = BENCHMARK_GROUPSIZE #hide
    params["gpu_groupsize_excitation"] = BENCHMARK_GROUPSIZE #hide
    return params #hide
end #hide

function benchmark_inputs(sim_method, gpu, receiver, ncoils) #hide
    obj = brain_phantom3D(; ss=3)[1:BENCHMARK_SPINS] #hide
    seq = PulseDesigner.EPI_example() #hide
    sys = Scanner(receiver=benchmark_receiver(receiver, ncoils)) #hide
    params = benchmark_parameters(sim_method, gpu) #hide
    return obj, seq, sys, params, gpu #hide
end #hide

function benchmark_simulation(obj, seq, sys, params, gpu) #hide
    signal = simulate(obj, seq, sys; sim_params=params, verbose=false) #hide
    benchmark_synchronize(gpu) #hide
    return signal #hide
end #hide

const BENCHMARK_CASES = ( #hide
    ( #hide
        key=:blochsimple_cpu, #hide
        label="BlochSimple CPU", #hide
        branch_uniform=2568.388, #hide
        master_uniform=2231.805, #hide
        birdcage=(2776.556, 3673.712, 4791.538), #hide
        arbitrary=(3368.451, 5038.578, 7091.253), #hide
    ), #hide
    ( #hide
        key=:bloch_cpu, #hide
        label="Bloch CPU", #hide
        branch_uniform=1070.796, #hide
        master_uniform=965.302, #hide
        birdcage=(2361.772, 6511.194, 12722.464), #hide
        arbitrary=(3755.384, 12992.677, 25834.163), #hide
    ), #hide
    ( #hide
        key=:bloch_gpu, #hide
        label="Bloch GPU", #hide
        branch_uniform=87.000, #hide
        master_uniform=68.885, #hide
        birdcage=(102.516, 183.089, 319.183), #hide
        arbitrary=(85.428, 190.842, 332.858), #hide
    ), #hide
    ( #hide
        key=:magnus1_cpu, #hide
        label="BlochMagnus1 CPU", #hide
        branch_uniform=842.647, #hide
        master_uniform=913.933, #hide
        birdcage=(2062.657, 6320.849, 12976.041), #hide
        arbitrary=(3884.993, 12853.516, 25857.504), #hide
    ), #hide
    ( #hide
        key=:magnus1_gpu, #hide
        label="BlochMagnus1 GPU", #hide
        branch_uniform=74.404, #hide
        master_uniform=67.576, #hide
        birdcage=(83.793, 183.696, 319.922), #hide
        arbitrary=(86.206, 192.619, 340.237), #hide
    ), #hide
    ( #hide
        key=:magnus2_cpu, #hide
        label="BlochMagnus2 CPU", #hide
        branch_uniform=978.603, #hide
        master_uniform=967.552, #hide
        birdcage=(2168.501, 6390.682, 12924.676), #hide
        arbitrary=(3800.820, 13010.392, 26136.354), #hide
    ), #hide
    ( #hide
        key=:magnus2_gpu, #hide
        label="BlochMagnus2 GPU", #hide
        branch_uniform=75.193, #hide
        master_uniform=68.261, #hide
        birdcage=(83.648, 183.943, 321.044), #hide
        arbitrary=(86.214, 191.517, 334.206), #hide
    ), #hide
    ( #hide
        key=:magnus4_cpu, #hide
        label="BlochMagnus4 CPU", #hide
        branch_uniform=956.589, #hide
        master_uniform=1042.492, #hide
        birdcage=(2126.222, 6280.580, 12886.550), #hide
        arbitrary=(3851.549, 12768.866, 25796.930), #hide
    ), #hide
    ( #hide
        key=:magnus4_gpu, #hide
        label="BlochMagnus4 GPU", #hide
        branch_uniform=74.569, #hide
        master_uniform=68.471, #hide
        birdcage=(85.057, 184.804, 318.806), #hide
        arbitrary=(86.159, 191.143, 333.521), #hide
    ), #hide
    ( #hide
        key=:magnus6_cpu, #hide
        label="BlochMagnus6 CPU", #hide
        branch_uniform=991.358, #hide
        master_uniform=1004.470, #hide
        birdcage=(2188.527, 6634.542, 14285.199), #hide
        arbitrary=(3896.261, 13470.822, 27440.335), #hide
    ), #hide
    ( #hide
        key=:magnus6_gpu, #hide
        label="BlochMagnus6 GPU", #hide
        branch_uniform=76.338, #hide
        master_uniform=76.365, #hide
        birdcage=(84.788, 183.348, 322.462), #hide
        arbitrary=(87.064, 191.628, 342.351), #hide
    ), #hide
) #hide

const BLOCHDICT_BENCHMARK_CASES = ( #hide
    ( #hide
        key=:blochdict_cpu, #hide
        label="BlochDict CPU", #hide
        uniform=18.981, #hide
        upstream_uniform=25.179, #hide
        birdcage=(16.716, 23.417, 27.492), #hide
        arbitrary=(19.865, 28.219, 41.157), #hide
    ), #hide
    ( #hide
        key=:blochdict_gpu, #hide
        label="BlochDict GPU", #hide
        uniform=25.448, #hide
        upstream_uniform=47.825, #hide
        birdcage=(19.550, 24.580, 31.607), #hide
        arbitrary=(24.224, 32.287, 45.421), #hide
    ), #hide
) #hide

function benchmark_plot(case, receiver) #hide
    receiver_label = String(receiver) #hide
    traces = [ #hide
        scatter( #hide
            x=[0, BENCHMARK_COILS...], #hide
            y=[case.branch_uniform, getproperty(case, receiver)...], #hide
            mode="lines+markers", #hide
            name="branch $(receiver_label)", #hide
            line=attr(color="#1f77b4", width=3), #hide
            marker=attr(size=9), #hide
            hovertemplate="coils=%{x}<br>time=%{y:.3f} ms<extra></extra>", #hide
        ), #hide
        scatter( #hide
            x=[0], #hide
            y=[case.master_uniform], #hide
            mode="markers", #hide
            name="master uniform", #hide
            marker=attr( #hide
                color="#111111", #hide
                size=14, #hide
                symbol="diamond-open", #hide
                line=attr(width=3), #hide
            ), #hide
            hovertemplate="master uniform<br>time=%{y:.3f} ms<extra></extra>", #hide
        ), #hide
    ] #hide
    return plot( #hide
        traces, #hide
        Layout( #hide
            title="$(case.label) with $(receiver_label) coils (100k spins)", #hide
            template="plotly_white", #hide
            width=950, #hide
            height=500, #hide
            margin=attr(l=60, r=30, t=70, b=60), #hide
            legend=attr(orientation="h", y=1.14, x=0.0), #hide
            xaxis=attr( #hide
                title="coil number", #hide
                tickmode="array", #hide
                tickvals=[0, BENCHMARK_COILS...], #hide
                ticktext=["uniform", string.(BENCHMARK_COILS)...], #hide
                range=[-0.25, 8.5], #hide
            ), #hide
            yaxis=attr(title="median time (ms)"), #hide
        ), #hide
        config=PlotConfig(scrollZoom=true), #hide
    ) #hide
end #hide

function blochdict_benchmark_plot(case, receiver) #hide
    receiver_label = String(receiver) #hide
    traces = [ #hide
        scatter( #hide
            x=[0, BENCHMARK_COILS...], #hide
            y=[case.uniform, getproperty(case, receiver)...], #hide
            mode="lines+markers", #hide
            name="branch $(receiver_label)", #hide
            line=attr(color="#1f77b4", width=3), #hide
            marker=attr(size=9), #hide
            hovertemplate="coils=%{x}<br>time=%{y:.3f} ms<extra></extra>", #hide
        ), #hide
        scatter( #hide
            x=[0], #hide
            y=[case.upstream_uniform], #hide
            mode="markers", #hide
            name="upstream uniform", #hide
            marker=attr( #hide
                color="#111111", #hide
                size=14, #hide
                symbol="diamond-open", #hide
                line=attr(width=3), #hide
            ), #hide
            hovertemplate="upstream uniform<br>time=%{y:.3f} ms<extra></extra>", #hide
        ), #hide
    ] #hide
    return plot( #hide
        traces, #hide
        Layout( #hide
            title="$(case.label) with $(receiver_label) coils (100k spins, 8 ADC)", #hide
            template="plotly_white", #hide
            width=950, #hide
            height=500, #hide
            margin=attr(l=60, r=30, t=70, b=60), #hide
            legend=attr(orientation="h", y=1.14, x=0.0), #hide
            xaxis=attr( #hide
                title="coil number", #hide
                tickmode="array", #hide
                tickvals=[0, BENCHMARK_COILS...], #hide
                ticktext=["uniform", string.(BENCHMARK_COILS)...], #hide
                range=[-0.25, 8.5], #hide
            ), #hide
            yaxis=attr(title="median time (ms)"), #hide
        ), #hide
        config=PlotConfig(scrollZoom=true), #hide
    ) #hide
end #hide

benchmark_plots = Dict( #hide
    (case.key, receiver) => benchmark_plot(case, receiver) #hide
    for case in BENCHMARK_CASES for receiver in (:birdcage, :arbitrary) #hide
) #hide

blochdict_benchmark_plots = Dict( #hide
    (case.key, receiver) => blochdict_benchmark_plot(case, receiver) #hide
    for case in BLOCHDICT_BENCHMARK_CASES for receiver in (:birdcage, :arbitrary) #hide
) #hide

function print_benchmark_summary() #hide
    for case in BENCHMARK_CASES #hide
        println("\n$(case.label)") #hide
        println("  branch uniform: $(case.branch_uniform) ms") #hide
        println("  master uniform: $(case.master_uniform) ms") #hide
        println( #hide
            "  branch uniform / master uniform: ", #hide
            round(case.branch_uniform / case.master_uniform; digits=3), #hide
            "x", #hide
        ) #hide
        for receiver in (:birdcage, :arbitrary), (index, ncoils) in pairs(BENCHMARK_COILS) #hide
            elapsed = getproperty(case, receiver)[index] #hide
            println( #hide
                "  $(receiver) $(ncoils) coil$(ncoils == 1 ? "" : "s"): ", #hide
                elapsed, #hide
                " ms  (", #hide
                round(elapsed / case.master_uniform; digits=3), #hide
                "x vs master uniform)", #hide
            ) #hide
        end #hide
    end #hide
    return nothing #hide
end; #hide

# ### BlochSimple CPU

# #### BirdcageCoilSens
#md benchmark_plots[(:blochsimple_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:blochsimple_cpu, :arbitrary)] #hide

# ### Bloch CPU

# #### BirdcageCoilSens
#md benchmark_plots[(:bloch_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:bloch_cpu, :arbitrary)] #hide

# ### Bloch GPU

# #### BirdcageCoilSens
#md benchmark_plots[(:bloch_gpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:bloch_gpu, :arbitrary)] #hide

# ### BlochMagnus1 CPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus1_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus1_cpu, :arbitrary)] #hide

# ### BlochMagnus1 GPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus1_gpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus1_gpu, :arbitrary)] #hide

# ### BlochMagnus2 CPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus2_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus2_cpu, :arbitrary)] #hide

# ### BlochMagnus2 GPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus2_gpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus2_gpu, :arbitrary)] #hide

# ### BlochMagnus4 CPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus4_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus4_cpu, :arbitrary)] #hide

# ### BlochMagnus4 GPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus4_gpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus4_gpu, :arbitrary)] #hide

# ### BlochMagnus6 CPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus6_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus6_cpu, :arbitrary)] #hide

# ### BlochMagnus6 GPU

# #### BirdcageCoilSens
#md benchmark_plots[(:magnus6_gpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md benchmark_plots[(:magnus6_gpu, :arbitrary)] #hide

# ## BlochDict Forward-Model Benchmark

# `BlochDict()` preserves every spin contribution instead of summing spins during
# acquisition. The 100,000-spin EPI benchmark above would therefore allocate a
# 7.6 GiB uniform dictionary and up to 60.8 GiB for eight receive-weighted
# channels. These BlochDict measurements retain the same 100,000-spin phantom,
# Float32 precision, receiver models, coil counts, CPU thread count, and Metal
# workgroup size, but use a separate RF excitation followed by eight ADC samples.
# The largest receive-weighted output is approximately 49 MiB.
#
# Each case measures dictionary simulation followed by static receive weighting.
# The values are medians of ten warmed `BenchmarkTools` samples with one
# evaluation per sample. CPU cases use 10 Julia threads; GPU cases use one CPU
# thread and 256 Metal work-items per workgroup. Their absolute times should not
# be compared with the full-EPI plots above.

# ### BlochDict CPU

# #### BirdcageCoilSens
#md blochdict_benchmark_plots[(:blochdict_cpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md blochdict_benchmark_plots[(:blochdict_cpu, :arbitrary)] #hide

# ### BlochDict GPU

# #### BirdcageCoilSens
#md blochdict_benchmark_plots[(:blochdict_gpu, :birdcage)] #hide

# #### ArbitraryCoilSens
#md blochdict_benchmark_plots[(:blochdict_gpu, :arbitrary)] #hide

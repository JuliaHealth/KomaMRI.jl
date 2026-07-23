# # Multi-Coil Receive Reconstruction

using KomaMRI, PlotlyJS, Suppressor #hide

# Receive coils measure the same transverse magnetization with different complex
# spatial weightings. This example begins with measured in-vitro sensitivities,
# then uses `Bloch()` on the CPU with arbitrary and birdcage receive models.

function coil_panel( #hide
    maps, label, colorscale, zmin, zmax; #hide
    columns=length(maps), show_axes=true, compact=false, #hide
) #hide
    ncoils = length(maps) #hide
    rows = cld(ncoils, columns) #hide
    titles = ["Coil $(coil) $(label)" for coil in 1:ncoils] #hide
    append!(titles, fill("", rows * columns - ncoils)) #hide
    panel = compact ? #hide
        make_subplots(; rows, cols=columns, horizontal_spacing=0.005, vertical_spacing=0.005) : #hide
        make_subplots(; #hide
            rows, cols=columns, #hide
            subplot_titles=reshape(titles, rows, columns), #hide
            horizontal_spacing=0.03, vertical_spacing=0.05, #hide
        ) #hide
    for (coil, map) in pairs(maps) #hide
        plot = plot_image(map; zmin, zmax, colorscale) #hide
        for trace in plot.data #hide
            trace.fields[:showscale] = false #hide
            add_trace!(panel, trace; row=cld(coil, columns), col=mod1(coil, columns)) #hide
        end #hide
        suffix = coil == 1 ? "" : string(coil) #hide
        xaxis = Symbol("xaxis", suffix) #hide
        yaxis = Symbol("yaxis", suffix) #hide
        xref = coil == 1 ? "x" : "x$(coil)" #hide
        Nx, Ny = size(map) #hide
        panel.plot.layout[xaxis][:range] = [-0.5, Nx - 0.5] #hide
        panel.plot.layout[yaxis][:range] = [-0.5, Ny - 0.5] #hide
        panel.plot.layout[yaxis][:scaleanchor] = xref #hide
        panel.plot.layout[yaxis][:constrain] = "domain" #hide
        panel.plot.layout[xaxis][:showticklabels] = show_axes #hide
        panel.plot.layout[yaxis][:showticklabels] = show_axes #hide
        panel.plot.layout[xaxis][:showgrid] = show_axes #hide
        panel.plot.layout[yaxis][:showgrid] = show_axes #hide
        panel.plot.layout[xaxis][:zeroline] = false #hide
        panel.plot.layout[yaxis][:zeroline] = false #hide
    end #hide
    if compact #hide
        relayout!( #hide
            panel; width="100%", height=135 * rows, showlegend=false, #hide
            margin=attr(l=0, r=0, t=0, b=0), #hide
            paper_bgcolor="#1b1b1f", plot_bgcolor="#000", #hide
        ) #hide
    else #hide
        relayout!(panel; width="100%", height=320 * rows, showlegend=false) #hide
    end #hide
    return panel.plot #hide
end; #hide

# ## Measured in-vitro coil sensitivities

# The scanner localizer is stored as ISMRMRD `.mrd` data. Its readout is
# oversampled, and its Cartesian metadata requires three corrections before
# reconstruction: restore the encoded phase dimension, match the encoded and
# reconstructed readout fields of view, and shift the phase-encode indices.

#md # ```julia
#md # using MAT, MRICoilSensitivities, MRIFiles, MRIReco
#md #
#md # raw = remove_oversampling(
#md #     RawAcquisitionData(ISMRMRDFile("localizer.mrd")),
#md # )
#md # raw.params["encodedSize"][2] = 256
#md # raw.params["encodedFOV"][1] = raw.params["reconFOV"][1]
#md # foreach(p -> p.head.idx.kspace_encode_step_1 += UInt16(12), raw.profiles)
#md # ```

# ESPIRiT uses a ``30 \times 30`` calibration region and a ``6 \times 6``
# kernel to estimate one complex sensitivity map for each of the 18 receive
# channels. The localizer contains three planes; the panels below show plane 1.

#md # ```julia
#md # maps = espirit(
#md #     AcquisitionData(raw), (6, 6), 30;
#md #     eigThresh_1=0.02, eigThresh_2=0.95,
#md # )
#md # matwrite(
#md #     "localizer_espirit_sensitivities.mat",
#md #     Dict("sensitivity_maps" => maps),
#md # )
#md # sensitivities_vitro = maps[:, :, 1, :]
#md # ```

# The complex ESPIRiT result is displayed directly. No additional KomaMRI
# forward simulation or EPI reconstruction is applied to these maps.

# Scanner channel `cha6` has no corresponding reference image, so its sixth grid
# position is left empty. This display omission does not remove the channel from
# ESPIRiT: all 18 estimated complex sensitivity maps are displayed below.

#md # ```@raw html
#md # <style>
#md # .vitro-grid{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:.35rem}
#md # .vitro-grid img{width:100%;aspect-ratio:1;object-fit:cover;object-position:center;display:block;background:#000}
#md # .vitro-empty{display:block;aspect-ratio:1;background:transparent}
#md # @media(max-width:900px){.vitro-grid{grid-template-columns:repeat(3,minmax(0,1fr))}}
#md # </style>
#md # <p><strong>Measured scanner reconstructions</strong></p>
#md # <div class="vitro-grid">
#md # <img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha0_slc0_ave0_ma.png" alt="Measured channel 0"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha1_slc0_ave0_ma.png" alt="Measured channel 1"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha2_slc0_ave0_ma.png" alt="Measured channel 2"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha3_slc0_ave0_ma.png" alt="Measured channel 3"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha4_slc0_ave0_ma.png" alt="Measured channel 4">
#md # <span class="vitro-empty" aria-label="Scanner channel 6 unavailable"></span>
#md # <img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha5_slc0_ave0_ma.png" alt="Measured channel 5"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha7_slc0_ave0_ma.png" alt="Measured channel 7"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha8_slc0_ave0_ma.png" alt="Measured channel 8"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha9_slc0_ave0_ma.png" alt="Measured channel 9"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha10_slc0_ave0_ma.png" alt="Measured channel 10"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha11_slc0_ave0_ma.png" alt="Measured channel 11">
#md # <img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha12_slc0_ave0_ma.png" alt="Measured channel 12"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha13_slc0_ave0_ma.png" alt="Measured channel 13"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha14_slc0_ave0_ma.png" alt="Measured channel 14"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha15_slc0_ave0_ma.png" alt="Measured channel 15"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha16_slc0_ave0_ma.png" alt="Measured channel 16"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha17_slc0_ave0_ma.png" alt="Measured channel 17">
#md # </div>
#md # ```

# ### ESPIRiT sensitivity magnitude

# Each map is independently normalized to emphasize its spatial receive profile.
# The 3-by-6 panel is interactive: hover to inspect values and drag to zoom.

vitro_shape = (256, 256, 18); #hide
vitro_data = joinpath( #hide
    dirname(pathof(KomaMRI)), "../docs/src/public/assets/", #hide
    "lit-10-vitro-coil-sens/espirit-magnitude-2d.f32", #hide
); #hide
vitro_reconstructions = open(vitro_data) do io #hide
    read!(io, Array{Float32}(undef, vitro_shape)) #hide
end; #hide
vitro_maps = [ #hide
    rot180(rotr90(permutedims(vitro_reconstructions[:, :, coil]))) #hide
    for coil in axes(vitro_reconstructions, 3) #hide
]; #hide
p_vitro = coil_panel( #hide
    vitro_maps, "magnitude", "Greys", 0, 1; #hide
    columns=6, show_axes=false, compact=true, #hide
) #hide
#md p_vitro #hide

# ### ESPIRiT sensitivity phase

# The phase of each complex ESPIRiT map is displayed inside its estimated object
# support. Colors span ``-\pi`` to ``\pi`` radians for all 18 channels.

vitro_phase_data = joinpath( #hide
    dirname(pathof(KomaMRI)), "../docs/src/public/assets/", #hide
    "lit-10-vitro-coil-sens/espirit-phase-2d.f32", #hide
); #hide
vitro_phase = open(vitro_phase_data) do io #hide
    read!(io, Array{Float32}(undef, vitro_shape)) #hide
end; #hide
vitro_phase_maps = [ #hide
    rot180(rotr90(permutedims(vitro_phase[:, :, coil]))) #hide
    for coil in axes(vitro_phase, 3) #hide
]; #hide
p_vitro_phase = coil_panel( #hide
    vitro_phase_maps, "phase", "Jet", -π, π; #hide
    columns=6, show_axes=false, compact=true, #hide
) #hide
#md p_vitro_phase #hide

# ## Synthetic acquisition setup

# The remaining examples use a 2D brain phantom and the same EPI sequence for
# both synthetic receiver models.

obj = brain_phantom2D();
seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq");
seq = @suppress read_seq(seq_file);

# The simulation parameters select the CPU `Bloch` method.

sim_params = KomaMRICore.default_sim_params();
sim_params["gpu"] = false;
sim_params["sim_method"] = Bloch();

# ## Arbitrary coil sensitivities

# `ArbitraryCoilSens` accepts sampled complex maps. Here, four Gaussian profiles
# model localized receive strength, while the complex phase varies around each
# coil center.

# In each equation, `exp` defines the Gaussian sensitivity magnitude. `complex`
# represents the direction from a coil center to the sample, `angle` extracts its
# angle, and `cis(θ) = cos(θ) + i sin(θ)` converts that angle into a unit-magnitude
# complex phase factor.

coil_spread = 0.08f0;
coil_offset = 0.20f0;

coil1(x, y, z) = exp(
    -Float32(π) * ((x - coil_offset)^2 + y^2 + z^2) / coil_spread,
) * cis(angle(complex(-y, -(x - coil_offset))));

coil2(x, y, z) = exp(
    -Float32(π) * (x^2 + (y - coil_offset)^2 + z^2) / coil_spread,
) * cis(angle(complex(-(y - coil_offset), -x)));

coil3(x, y, z) = exp(
    -Float32(π) * ((x + coil_offset)^2 + y^2 + z^2) / coil_spread,
) * cis(angle(complex(-y, -(x + coil_offset))));

coil4(x, y, z) = exp(
    -Float32(π) * (x^2 + (y + coil_offset)^2 + z^2) / coil_spread,
) * cis(angle(complex(-(y + coil_offset), -x)));

# We sample those profiles on a 3D grid and store the coil index in the fourth
# dimension of `coil_sens`.

coords = LinRange(-0.23f0, 0.23f0, 17);
zcoords = LinRange(-0.01f0, 0.01f0, 3);
coil_sens1 = [coil1(x, y, z) for x in coords, y in coords, z in zcoords];
coil_sens2 = [coil2(x, y, z) for x in coords, y in coords, z in zcoords];
coil_sens3 = [coil3(x, y, z) for x in coords, y in coords, z in zcoords];
coil_sens4 = [coil4(x, y, z) for x in coords, y in coords, z in zcoords];
coil_sens = cat(coil_sens1, coil_sens2, coil_sens3, coil_sens4; dims=4);
arbitrary = ArbitraryCoilSens(coords, coords, zcoords, coil_sens);

# The acquisition steps are unchanged; only the receiver attached to the scanner
# differs.

sys_arbitrary = Scanner(receiver=arbitrary);
raw_arbitrary = @suppress simulate(
    obj, seq, sys_arbitrary; sim_params, verbose=false
);

# We reconstruct the four arbitrary-coil signals in the same way.

acq_arbitrary = AcquisitionData(raw_arbitrary);
acq_arbitrary.traj[1].circular = false; #hide
Nx, Ny = raw_arbitrary.params["reconSize"][1:2];
recon_params = Dict(:reco => "direct", :reconSize => (Nx, Ny));
image_arbitrary = @suppress reconstruction(acq_arbitrary, recon_params);

# The physical reconstruction coordinates are used below to evaluate the
# sensitivity maps at every image pixel.

FOVx, FOVy = raw_arbitrary.params["reconFOV"][1:2] .* 1f-3;
xs = range(-FOVx / 2, FOVx / 2; length=Nx);
ys = range(-FOVy / 2, FOVy / 2; length=Ny);
xgrid = vec([x for x in xs, _ in ys]);
ygrid = vec([y for _ in xs, y in ys]);
zgrid = zeros(Float32, length(xgrid));

# ### Coil reconstructions

# The arbitrary-coil reconstruction has the same dimensions. We again retain all
# image pixels and coil channels, then extract one magnitude image per coil.

arbitrary_images = image_arbitrary[:, :, 1, 1, :, 1];
arbitrary_recon = [
    abs.(arbitrary_images[:, :, coil]) for coil in axes(arbitrary_images, 3)
];

p4 = coil_panel( #hide
    arbitrary_recon, "reco", "Greys", #hide
    minimum(minimum, arbitrary_recon), maximum(maximum, arbitrary_recon), #hide
) #hide
#md p4 #hide

# ### Sensitivity magnitude

# `get_sens` now interpolates the sampled 3D maps onto the same reconstruction
# coordinates. It again returns one complex column per coil. Taking `abs.` and
# reshaping each column produces the four sensitivity-magnitude images.

arbitrary_sens = get_sens(arbitrary, xgrid, ygrid, zgrid);
arbitrary_magnitude = [
    reshape(abs.(arbitrary_sens[:, coil]), Nx, Ny)
    for coil in axes(arbitrary_sens, 2)
];

p5 = coil_panel( #hide
    arbitrary_magnitude, "magnitude", "Greys", #hide
    minimum(minimum, arbitrary_magnitude), maximum(maximum, arbitrary_magnitude), #hide
) #hide
#md p5 #hide

# ### Sensitivity phase

# Finally, `angle.` extracts the interpolated complex phase in radians and
# `reshape` maps each flattened coil column back onto the image grid.

arbitrary_phase = [
    reshape(angle.(arbitrary_sens[:, coil]), Nx, Ny)
    for coil in axes(arbitrary_sens, 2)
];

p6 = coil_panel(arbitrary_phase, "phase", "Jet", -π, π) #hide
#md p6 #hide

# ## Birdcage coil sensitivities

# `BirdcageCoilSens` places four idealized receive elements uniformly around the
# object. For channel ``c``, the received signal is weighted by its complex
# sensitivity ``S_c(\mathbf{r})``. The number of coils, coil-ring radius, and
# axial half-length can be changed with `ncoils`, `radius`, and `L`, respectively.

birdcage = BirdcageCoilSens(ncoils=4, radius=0.18, L=0.20);
sys_birdcage = Scanner(receiver=birdcage);
raw_birdcage = @suppress simulate(
    obj, seq, sys_birdcage; sim_params, verbose=false);

# We convert the raw signals to acquisition data and reconstruct every receive
# channel directly.

acq_birdcage = AcquisitionData(raw_birdcage);
acq_birdcage.traj[1].circular = false; #hide
Nx, Ny = raw_birdcage.params["reconSize"][1:2];
recon_params = Dict(:reco => "direct", :reconSize => (Nx, Ny));
image_birdcage = @suppress reconstruction(acq_birdcage, recon_params);

# ### Coil reconstructions

# The reconstruction has dimensions
# `x × y × slice × contrast × coil × repetition`. This acquisition contains one
# slice, one contrast, and one repetition, so we select index `1` along those
# dimensions while retaining every image pixel and all four coil channels. The
# result is an `Nx × Ny × 4` array of complex-valued coil images.

birdcage_images = image_birdcage[:, :, 1, 1, :, 1];

# Each complex coil image contains magnitude and phase. Here we extract the
# magnitude with `abs.` and store the four resulting `Nx × Ny` images separately
# for plotting.

birdcage_recon = [
    abs.(birdcage_images[:, :, coil]) for coil in axes(birdcage_images, 3)
];

p1 = coil_panel( #hide
    birdcage_recon, "reco", "Greys", #hide
    minimum(minimum, birdcage_recon), maximum(maximum, birdcage_recon), #hide
) #hide
#md p1 #hide

# ### Sensitivity magnitude

# These maps describe the spatial receive strength of each channel. We first
# construct the physical location of every reconstruction pixel from the field of
# view and image dimensions. `xgrid`, `ygrid`, and `zgrid` are flattened because
# `get_sens` expects one ``(x,y,z)`` coordinate per pixel.

Nx, Ny = raw_birdcage.params["reconSize"][1:2];
FOVx, FOVy = raw_birdcage.params["reconFOV"][1:2] .* 1f-3;
xs = range(-FOVx / 2, FOVx / 2; length=Nx);
ys = range(-FOVy / 2, FOVy / 2; length=Ny);
xgrid = vec([x for x in xs, _ in ys]);
ygrid = vec([y for _ in xs, y in ys]);
zgrid = zeros(Float32, length(xgrid));

# `get_sens` evaluates the birdcage model at those coordinates and returns a
# complex `(Nx * Ny) × ncoils` matrix. Each column contains one coil's
# sensitivities. `abs.` extracts their magnitudes, and `reshape` restores each
# flattened column to an `Nx × Ny` image. `axes(..., 2)` iterates over the coil
# columns.

birdcage_sens = get_sens(birdcage, xgrid, ygrid, zgrid);
birdcage_magnitude = [
    reshape(abs.(birdcage_sens[:, coil]), Nx, Ny)
    for coil in axes(birdcage_sens, 2)
];

p2 = coil_panel( #hide
    birdcage_magnitude, "magnitude", "Greys", #hide
    minimum(minimum, birdcage_magnitude), maximum(maximum, birdcage_magnitude), #hide
) #hide
#md p2 #hide

# ### Sensitivity phase

# These maps show the spatial phase contributed by each channel. `angle.`
# extracts the phase of every complex sensitivity in radians, and `reshape`
# restores each coil's phase samples to the reconstruction grid. The plots use a
# common range from ``-\pi`` to ``\pi`` so colors represent the same phase in all
# four channels.

birdcage_phase = [
    reshape(angle.(birdcage_sens[:, coil]), Nx, Ny)
    for coil in axes(birdcage_sens, 2)
];

p3 = coil_panel(birdcage_phase, "phase", "Jet", -π, π) #hide
#md p3 #hide

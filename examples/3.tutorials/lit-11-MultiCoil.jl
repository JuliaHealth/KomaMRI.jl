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
        for trace in plot.plot.data #hide
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
    return panel #hide
end; #hide

# ## Measured in-vitro coil sensitivities

# A scanner calibration acquisition can provide measured receive sensitivities.
# Vendor raw data, such as a Siemens `.dat` file, is first converted to ISMRMRD
# `.mrd`. `MRIFiles.jl` loads the data, `MRIReco.jl` reconstructs each receive
# channel, and `MRICoilSensitivities.jl` estimates one complex 3D sensitivity map
# per channel. For this fully sampled calibration scan, the maps are obtained by
# normalizing the complex coil images by their root-sum-of-squares image.

#md # ```julia
#md # using MRIReco, MRIFiles, MRICoilSensitivities
#md #
#md # raw_vitro = RawAcquisitionData(ISMRMRDFile("coil_calibration.mrd"));
#md # ncoils_vitro = raw_vitro.params["receiverChannels"];
#md # profiles_vitro = filter(raw_vitro.profiles) do profile
#md #     size(profile.data, 2) == ncoils_vitro &&
#md #         flag_is_set(profile, "ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA")
#md # end;
#md # acq_vitro = AcquisitionData(RawAcquisitionData(raw_vitro.params, profiles_vitro));
#md # encoded_size = Tuple(raw_vitro.params["encodedSize"]);
#md # image_vitro = reconstruction(
#md #     acq_vitro, Dict(:reco => "direct", :reconSize => encoded_size),
#md # );
#md # sensitivities_vitro = estimateCoilSensitivities(image_vitro);
#md # ```

# `espirit(acq_vitro, (4, 4, 4), 16)` provides an alternative for
# autocalibration data. Here, readout oversampling is removed from the direct
# estimate before reordering the MRD axes into physical ``(x,y,z)`` order,
# attaching physical coordinates, and constructing `ArbitraryCoilSens`.

#md # ```julia
#md # nread = raw_vitro.params["reconSize"][1];
#md # first_read = (size(sensitivities_vitro, 1) - nread) ÷ 2 + 1;
#md # read_range = first_read:(first_read + nread - 1);
#md # coil_maps_vitro = sensitivities_vitro[read_range, :, :, 1, :, 1];
#md #
#md # fov_vitro = Float32.(raw_vitro.params["reconFOV"]) .* 1f-3;
#md # x_vitro = range(-fov_vitro[1] / 2, fov_vitro[1] / 2; length=size(coil_maps_vitro, 1));
#md # y_vitro = range(-fov_vitro[2] / 2, fov_vitro[2] / 2; length=size(coil_maps_vitro, 2));
#md # z_vitro = range(-fov_vitro[3] / 2, fov_vitro[3] / 2; length=size(coil_maps_vitro, 3));
#md # coil_magnitudes_vitro = complex.(abs.(coil_maps_vitro));
#md # receiver_vitro = ArbitraryCoilSens(
#md #     x_vitro, y_vitro, z_vitro, coil_magnitudes_vitro,
#md # );
#md # z0 = argmin(abs.(receiver_vitro.z));
#md # coil1_2D = receiver_vitro.coil_sens[:, :, z0, 1];
#md # coil2_2D = receiver_vitro.coil_sens[:, :, z0, 2];
#md # coil3_2D = receiver_vitro.coil_sens[:, :, z0, 3];
#md # coil4_2D = receiver_vitro.coil_sens[:, :, z0, 4];
#md # ```

# The in-vitro object is spherical, so this comparison models its central slice
# as a uniform disk rather than using the brain phantom shown later. The scanner
# reference images contain magnitude only, so this comparison uses the measured
# sensitivity magnitudes with zero phase. This preserves each coil's spatial
# receive weighting without introducing the unresolved phase of the estimated
# maps. Equal-area golden-angle positions approximate a continuous object, the
# sphere is simulated with `Bloch()`, and all 18 channels are reconstructed
# independently.

#md # ```julia
#md # radius = 50f-3;
#md # nspins = 100_000;
#md # spin_index = 1:nspins;
#md # spin_radius = @. radius * sqrt((spin_index - 0.5f0) / nspins);
#md # spin_angle = @. Float32(π) * (3 - sqrt(5f0)) * spin_index;
#md # spin_x = spin_radius .* cos.(spin_angle);
#md # spin_y = spin_radius .* sin.(spin_angle);
#md # obj_vitro = Phantom(;
#md #     name="uniform sphere slice",
#md #     x=spin_x,
#md #     y=spin_y,
#md # );
#md #
#md # seq_file_vitro = joinpath(
#md #     dirname(pathof(KomaMRI)),
#md #     "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq",
#md # );
#md # seq_vitro = read_seq(seq_file_vitro);
#md # sys_vitro = Scanner(receiver=receiver_vitro);
#md # sim_params_vitro = KomaMRICore.default_sim_params();
#md # sim_params_vitro["gpu"] = false;
#md # sim_params_vitro["sim_method"] = Bloch();
#md # raw_sim_vitro = simulate(
#md #     obj_vitro, seq_vitro, sys_vitro; sim_params=sim_params_vitro,
#md # );
#md # acq_sim_vitro = AcquisitionData(raw_sim_vitro);
#md # acq_sim_vitro.traj[1].circular = false;
#md # image_sim_vitro = reconstruction(
#md #     acq_sim_vitro, Dict(:reco => "direct", :reconSize => (100, 100)),
#md # );
#md # ```

# Both rows below use the central `z` slice and receiver channels `1:18` from
# the same MRD dataset. Each image is independently normalized so the comparison
# emphasizes spatial receive weighting rather than absolute signal intensity.

#md # ```@raw html
#md # <style>
#md # .vitro-grid{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:.35rem}
#md # .vitro-grid img{width:100%;aspect-ratio:1;object-fit:cover;object-position:center;display:block;background:#000}
#md # @media(max-width:900px){.vitro-grid{grid-template-columns:repeat(3,minmax(0,1fr))}}
#md # </style>
#md # <p><strong>Measured scanner reconstructions</strong></p>
#md # <div class="vitro-grid">
#md # <img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha0_slc0_ave0_ma.png" alt="Measured reconstruction 1"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha1_slc0_ave0_ma.png" alt="Measured reconstruction 2"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha2_slc0_ave0_ma.png" alt="Measured reconstruction 3"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha3_slc0_ave0_ma.png" alt="Measured reconstruction 4"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha4_slc0_ave0_ma.png" alt="Measured reconstruction 5"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha4_slc1_ave0_ma.png" alt="Measured reconstruction 6">
#md # <img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha5_slc0_ave0_ma.png" alt="Measured reconstruction 7"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha7_slc0_ave0_ma.png" alt="Measured reconstruction 8"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha8_slc0_ave0_ma.png" alt="Measured reconstruction 9"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha9_slc0_ave0_ma.png" alt="Measured reconstruction 10"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha10_slc0_ave0_ma.png" alt="Measured reconstruction 11"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha11_slc0_ave0_ma.png" alt="Measured reconstruction 12">
#md # <img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha12_slc0_ave0_ma.png" alt="Measured reconstruction 13"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha13_slc0_ave0_ma.png" alt="Measured reconstruction 14"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha14_slc0_ave0_ma.png" alt="Measured reconstruction 15"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha15_slc0_ave0_ma.png" alt="Measured reconstruction 16"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha16_slc0_ave0_ma.png" alt="Measured reconstruction 17"><img src="/assets/lit-10-vitro-coil-sens/MID0_COLxLIN_cha17_slc0_ave0_ma.png" alt="Measured reconstruction 18">
#md # </div>
#md # ```

# ### KomaMRI reconstructions using measured `ArbitraryCoilSens`

# The simulated images use the same physical display field of view as the
# measured center crop. Unlike the scanner PNGs above, this 3-by-6 panel is
# interactive: hover to inspect normalized magnitude and drag to zoom.

vitro_shape = (56, 56, 18); #hide
vitro_data = joinpath( #hide
    dirname(pathof(KomaMRI)), "../docs/src/public/assets/", #hide
    "lit-10-vitro-coil-sens/simulated-magnitude-axial.f32", #hide
); #hide
vitro_reconstructions = open(vitro_data) do io #hide
    read!(io, Array{Float32}(undef, vitro_shape)) #hide
end; #hide
vitro_maps = [ #hide
    permutedims(vitro_reconstructions[:, :, coil]) #hide
    for coil in axes(vitro_reconstructions, 3) #hide
]; #hide
p_vitro = coil_panel( #hide
    vitro_maps, "reco", "Greys", 0, 1; #hide
    columns=6, show_axes=false, compact=true, #hide
) #hide
#md p_vitro #hide

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

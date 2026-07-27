# # Parallel Imaging with Receive Coils

using KomaMRI, PlotlyBase, Suppressor #hide

# Parallel imaging acquires fewer phase-encoding lines and uses the distinct
# spatial sensitivities of multiple receive coils to recover the missing
# information. We will compare fully sampled and two-fold accelerated EPI
# acquisitions of the same 2D brain phantom.

fat_free_tissues = Dict(
    "FAT1" => zeros(5),
    "FAT2" => zeros(5),
);
obj = brain_phantom2D(; tissue_properties=fat_free_tissues);

# ## Fully sampled and accelerated EPI

# We begin with a 100-by-100 EPI sequence. The accelerated sequence keeps half
# of its readouts and doubles the phase-encoding blips, so it traverses every
# second ``k_y`` line with acceleration factor ``R=2``.

seq_file = joinpath(
    dirname(pathof(KomaMRI)),
    "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/",
    "epi_100x100_TE100_FOV230.seq",
);
full_seq = @suppress read_seq(seq_file);

R = 2;
readout_blocks = findall(i -> is_ADC_on(full_seq[i]), 1:length(full_seq));
last_accelerated_readout = readout_blocks[length(readout_blocks) ÷ R];
trailing_blocks = (last(readout_blocks) + 1):length(full_seq);
accelerated_seq = deepcopy(
    full_seq[1:last_accelerated_readout] + full_seq[trailing_blocks],
);
phase_blips = (first(readout_blocks) + 1):2:(last_accelerated_readout - 1);
accelerated_seq.GR[2, phase_blips] .= R .* accelerated_seq.GR[2, phase_blips];
_, full_kspace = get_kspace(full_seq);
_, accelerated_kspace = get_kspace(accelerated_seq);
full_center_sample = argmin(vec(sum(abs2, full_kspace; dims=2)));
accelerated_center_sample = argmin(vec(sum(abs2, accelerated_kspace; dims=2)));
echo_time_delay =
    get_adc_sampling_times(full_seq)[full_center_sample] -
    get_adc_sampling_times(accelerated_seq)[accelerated_center_sample];
accelerated_seq =
    accelerated_seq[1:(first(readout_blocks) - 1)] +
    Delay(echo_time_delay) +
    accelerated_seq[first(readout_blocks):end];
accelerated_seq.DEF["Name"] = "epi_R2";
accelerated_seq.DEF["Num_Blocks"] = length(accelerated_seq);

# `plot_kspace(...; view_2d=true)` makes the missing phase-encoding lines
# explicit. The left trajectory has 100 acquired lines; the right has 50. The
# added delay keeps the time from excitation to the centre of k-space equal in
# both sequences.

p_kspace_full = plot_kspace( #hide
    full_seq; view_2d=true, width=430, height=430, #hide
) #hide
p_kspace_accelerated = plot_kspace( #hide
    accelerated_seq; view_2d=true, width=430, height=430, #hide
) #hide
relayout!(p_kspace_full; title="Fully sampled k-space") #hide
relayout!(p_kspace_accelerated; title="R = 2 k-space") #hide
foreach(trace -> trace.fields[:showlegend] = false, p_kspace_accelerated.data) #hide

# ## Without parallel imaging

# First we simulate both sequences with the default single uniform receive
# channel. A direct reconstruction has no information with which to separate
# pixels that overlap after the ``R=2`` undersampling.

uniform_sys = Scanner();
raw_uniform_full = @suppress simulate(
    obj, full_seq, uniform_sys; verbose=false,
);
raw_uniform_accelerated = @suppress simulate(
    obj, accelerated_seq, uniform_sys; verbose=false,
);

acq_uniform_full = AcquisitionData(raw_uniform_full);
acq_uniform_accelerated = AcquisitionData(raw_uniform_accelerated);
acq_uniform_full.traj[1].circular = false;
acq_uniform_accelerated.traj[1].circular = false;

recon_size = Tuple(raw_uniform_full.params["reconSize"][1:2]);
direct_params = Dict{Symbol,Any}(
    :reco => "direct",
    :reconSize => recon_size,
);
image_uniform_full = reconstruction(acq_uniform_full, direct_params);
image_uniform_accelerated = reconstruction(
    acq_uniform_accelerated, direct_params,
);

p_uniform_full = plot_image( #hide
    abs.(image_uniform_full[:, :, 1, 1, 1, 1]); #hide
    width=430, height=430, title="Fully sampled", #hide
) #hide
p_uniform_accelerated = plot_image( #hide
    abs.(image_uniform_accelerated[:, :, 1, 1, 1, 1]); #hide
    width=430, height=430, title="R = 2: direct reconstruction", #hide
) #hide
p_uniform_full.data[1].fields[:showscale] = false #hide
p_uniform = [ #hide
    p_kspace_full p_kspace_accelerated; #hide
    p_uniform_full p_uniform_accelerated #hide
] #hide
relayout!(p_uniform; width=950, height=900) #hide
#md p_uniform #hide
#jl display(p_uniform)

# The accelerated reconstruction contains two overlapping copies of the brain,
# separated by half the field of view. This is the aliasing that the receive
# array must resolve.

# ## Add a birdcage receive array

# We repeat the simulations with an eight-channel birdcage receiver. The pulse
# sequence and phantom are unchanged; only the receiver attached to the scanner
# differs.

receiver = BirdcageCoilSens(; ncoils=8);
coil_sys = Scanner(; receiver);
raw_coils_full = @suppress simulate(
    obj, full_seq, coil_sys; verbose=false,
);
raw_coils_accelerated = @suppress simulate(
    obj, accelerated_seq, coil_sys; verbose=false,
);

acq_coils_full = AcquisitionData(raw_coils_full);
acq_coils_accelerated = AcquisitionData(raw_coils_accelerated);
acq_coils_full.traj[1].circular = false;
acq_coils_accelerated.traj[1].circular = false;

# SENSE needs one complex sensitivity value for every reconstruction pixel and
# channel. We therefore evaluate the receiver on the physical reconstruction
# grid, not on the phantom's irregular spin locations.

Nx, Ny = recon_size;
FOVx, FOVy = raw_coils_full.params["reconFOV"][1:2] .* 1e-3;
T = typeof(real(zero(eltype(first(raw_coils_full.profiles).data))));
x_axis = range(T(-FOVx / 2), T(FOVx / 2); length=Nx);
y_axis = range(T(-FOVy / 2), T(FOVy / 2); length=Ny);
x_positions = vec([x for x in x_axis, _ in y_axis]);
y_positions = vec([y for _ in x_axis, y in y_axis]);
z_positions = zeros(T, length(x_positions));
sensitivity_maps = reshape(
    get_sens(receiver, x_positions, y_positions, z_positions),
    Nx, Ny, 1, get_n_coils(receiver),
);

# `BirdcageCoilSens` evaluates an analytic model. For measured maps, an
# `ArbitraryCoilSens` receiver uses the same `get_sens` call and linearly
# interpolates its stored Cartesian samples onto this reconstruction grid.

# A direct reconstruction of the fully sampled multi-channel data retains one
# image per receive channel. The first row below shows the magnitude of the
# first four sensitivity maps; the second row shows their corresponding coil
# images. Each sensitivity is normalized independently to emphasize its spatial
# profile.

coil_images = reconstruction(acq_coils_full, direct_params);
displayed_coils = 1:4;
coil_image_scale = maximum(abs, coil_images[:, :, 1, 1, displayed_coils, 1]); #hide
sensitivity_plots = [ #hide
    plot_image( #hide
        abs.(sensitivity_maps[:, :, 1, coil]) ./ #hide
        maximum(abs, sensitivity_maps[:, :, 1, coil]); #hide
        width=220, height=260, zmin=0, zmax=1, #hide
        title="Coil $coil sensitivity", #hide
    ) #hide
    for coil in displayed_coils #hide
] #hide
coil_image_plots = [ #hide
    plot_image( #hide
        abs.(coil_images[:, :, 1, 1, coil, 1]); #hide
        width=220, height=260, zmin=0, zmax=coil_image_scale, #hide
        title="Coil $coil image", #hide
    ) #hide
    for coil in displayed_coils #hide
] #hide
foreach( #hide
    plot -> plot.data[1].fields[:showscale] = false, #hide
    [sensitivity_plots; coil_image_plots], #hide
) #hide
p_coils = [ #hide
    sensitivity_plots[1] sensitivity_plots[2] sensitivity_plots[3] sensitivity_plots[4]; #hide
    coil_image_plots[1] coil_image_plots[2] coil_image_plots[3] coil_image_plots[4] #hide
] #hide
relayout!(p_coils; width=1100, height=600) #hide
#md p_coils #hide
#jl display(p_coils)

# ## SENSE reconstruction

# The SENSE forward model combines the image with these complex maps before
# sampling its k-space trajectory. MRIReco's multi-coil reconstruction solves
# the inverse problem jointly across all eight channels.

sense_params = Dict{Symbol,Any}(
    :reco => "multiCoil",
    :reconSize => recon_size,
    :senseMaps => sensitivity_maps,
    :iterations => 20,
    :densityWeighting => false,
    :toeplitz => false,
);
sense_full = reconstruction(acq_coils_full, sense_params);
sense_accelerated = reconstruction(acq_coils_accelerated, sense_params);

p_sense_full = plot_image( #hide
    abs.(sense_full[:, :, 1, 1, 1, 1]); #hide
    width=430, height=430, title="Fully sampled SENSE", #hide
) #hide
p_sense_accelerated = plot_image( #hide
    abs.(sense_accelerated[:, :, 1, 1, 1, 1]); #hide
    width=430, height=430, title="R = 2 SENSE", #hide
) #hide
p_sense_full.data[1].fields[:showscale] = false #hide
p_sense = [p_sense_full p_sense_accelerated] #hide
relayout!(p_sense; width=950, height=450) #hide
#md p_sense #hide
#jl display(p_sense)

# The accelerated SENSE image no longer contains the fold-over seen with the
# single receive channel. The missing phase-encoding lines are recovered from
# the spatially distinct complex measurements supplied by the receiver array.

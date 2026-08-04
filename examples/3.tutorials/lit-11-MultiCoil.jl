# # Accelerated MRI with SENSE

using KomaMRI, PlotlyBase #hide
#hide
function square_axes!(plot; hide_ticks=false, height) #hide
    for (name, axis) in plot.layout.fields #hide
        axis_name = String(name) #hide
        if startswith(axis_name, "yaxis") #hide
            suffix = replace(axis_name, "yaxis" => "") #hide
            axis[:scaleanchor] = suffix in ("", "1") ? "x" : "x$suffix" #hide
            axis[:scaleratio] = 1 #hide
            axis[:constrain] = "domain" #hide
        end #hide
        if hide_ticks && (startswith(axis_name, "xaxis") || startswith(axis_name, "yaxis")) #hide
            axis[:showticklabels] = false #hide
            axis[:title] = attr(text="") #hide
            axis[:showgrid] = false #hide
            axis[:zeroline] = false #hide
        end #hide
    end #hide
    relayout!(plot; height, margin=attr(l=20, r=20, t=55, b=20)) #hide
    return plot #hide
end #hide
#hide
function image_panel(image, title; zmax) #hide
    plot = plot_image(image; title, zmin=0, zmax) #hide
    plot.data[1].fields[:showscale] = false #hide
    return plot #hide
end #hide
#hide
function direct_reconstruction(raw, recon_size) #hide
    acquisition = AcquisitionData(raw) #hide
    acquisition.traj[1].circular = false #hide
    params = Dict{Symbol,Any}(:reco => "direct", :reconSize => recon_size) #hide
    return reconstruction(acquisition, params), acquisition #hide
end #hide
#hide
tissue_properties = Dict("FAT1" => zeros(5), "FAT2" => zeros(5)); #hide
obj = brain_phantom2D(; tissue_properties); #hide
nothing #hide

# Parallel imaging reduces scan time by acquiring fewer phase-encoding lines.
# The missing data make a conventional image fold over, but an array of receive
# coils provides additional spatial information that can separate the
# overlapping pixels. We will compare fully sampled and three-fold accelerated
# EPI acquisitions of the same
# [`brain_phantom2D`](@ref KomaMRIBase.brain_phantom2D).

# ## Undersample k-space

# The fully sampled EPI trajectory visits every ``k_y`` line. The accelerated
# sequence visits every third line, so its acceleration factor is ``R=3``.
# Both sequences reach the centre of k-space at the same echo time. The
# [`plot_kspace`](@ref KomaMRIPlots.plot_kspace) view below makes the missing
# lines explicit.

function epi_sequence(R; matrix_size=120, FOV=0.23, TE=50e-3, dwell=4e-6, sys=Scanner()) #hide
    readout = PulseDesigner.make_trapezoid( #hide
        ; flat_area=matrix_size / FOV, flat_time=matrix_size * dwell, sys, #hide
    ) #hide
    adc = PulseDesigner.make_adc( #hide
        matrix_size; dwell, delay=readout.rise, sys, #hide
    ) #hide
    ky_lines = (-matrix_size ÷ 2):R:(matrix_size ÷ 2 - 1) #hide
    train = Sequence(sys) #hide
    @addblock train += ( #hide
        x=PulseDesigner.make_trapezoid(; area=-γ * area(readout) / 2, sys), #hide
        y=PulseDesigner.make_trapezoid(; area=first(ky_lines) / FOV, sys), #hide
    ) #hide
    blip = PulseDesigner.make_trapezoid(; area=R / FOV, sys) #hide
    for (line, _) in enumerate(ky_lines) #hide
        @addblock train += (x=isodd(line) ? readout : -readout, adc) #hide
        if line < length(ky_lines) #hide
            @addblock train += (y=blip) #hide
        end #hide
    end #hide
    center_block = 2 * findfirst(iszero, ky_lines) #hide
    center_time = get_block_start_times(train)[center_block] + #hide
        adc.delay + adc.T / 2 #hide
    excitation = PulseDesigner.build_block_pulse( #hide
        π / 2; duration=1e-3, use=Excitation(), sys, #hide
    ) #hide
    rf = excitation.RF[1, 1] #hide
    excitation_center = rf.delay + rf_center(rf) #hide
    pre_readout_delay = PulseDesigner.round_to_raster( #hide
        TE - (dur(excitation) - excitation_center + center_time), #hide
        sys.limits.DUR_Δt, #hide
    ) #hide
    seq = Sequence(sys) #hide
    @addblock seq += excitation #hide
    @addblock seq += Delay(pre_readout_delay) #hide
    @addblock seq += train #hide
    return seq #hide
end #hide
#hide
R = 3 #hide
matrix_size = 120 #hide
FOV = 0.23 #hide
full_seq = epi_sequence(1; matrix_size, FOV) #hide
acc_seq = epi_sequence(R; matrix_size, FOV) #hide
_, full_kspace = get_kspace(full_seq) #hide
_, acc_kspace = get_kspace(acc_seq) #hide
#hide
p_kspace_full = plot_kspace(full_seq; view_2d=true) #hide
p_kspace_acc = plot_kspace( #hide
    acc_seq; view_2d=true, #hide
) #hide
relayout!(p_kspace_full; title="Fully sampled") #hide
relayout!(p_kspace_acc; title="R = 3") #hide
foreach( #hide
    trace -> trace.fields[:showlegend] = false, #hide
    [p_kspace_full.data; p_kspace_acc.data], #hide
) #hide
p_kspace = square_axes!([p_kspace_full p_kspace_acc]; height=300) #hide
#md p_kspace #hide
#jl display(p_kspace)

# ## See the aliasing with sum-of-squares

# Attach a 16-channel
# [`BirdcageCoilSens`](@ref KomaMRIBase.BirdcageCoilSens) receiver. The
# sequence and phantom do not change; each channel simply measures the same
# transverse magnetization through a different complex spatial sensitivity.

birdcage = BirdcageCoilSens(; ncoils=16);
system = Scanner(; receiver=birdcage);
raw_coils_full = simulate(obj, full_seq, system; verbose=false);
raw_coils_acc = simulate(obj, acc_seq, system; verbose=false);
#hide
recon_size = (matrix_size, matrix_size) #hide
coil_images_full, _ = direct_reconstruction(raw_coils_full, recon_size); #hide
coil_images_acc, acq_coils_acc = direct_reconstruction(raw_coils_acc, recon_size); #hide
#hide
# Sum-of-squares (SoS) combines the reconstructed channel magnitudes. It
# produces a robust fully sampled image, but it cannot separate pixels that
# overlap after undersampling:

sos_full = sqrt.(sum(abs2, coil_images_full; dims=5));
sos_acc = sqrt.(sum(abs2, coil_images_acc; dims=5));

sos_full_image = sos_full[:, :, 1, 1, 1, 1] #hide
sos_acc_image = sos_acc[:, :, 1, 1, 1, 1] #hide
p_sos_full = image_panel( #hide
    sos_full_image, "Fully sampled SoS"; zmax=maximum(sos_full_image), #hide
) #hide
p_sos_acc = image_panel( #hide
    sos_acc_image, "R = 3 SoS: aliased"; zmax=maximum(sos_acc_image), #hide
) #hide
p_sos = square_axes!( #hide
    [p_sos_full p_sos_acc]; #hide
    hide_ticks=true, height=300, #hide
) #hide
#md p_sos #hide
#jl display(p_sos)

# ## Inspect the coil sensitivities

# SENSE needs one complex sensitivity per reconstruction pixel and receive
# channel. Use [`get_sens`](@ref KomaMRIBase.get_sens) to evaluate the receiver
# on the reconstruction grid—not on the phantom's irregular spin positions—and
# [`get_n_coils`](@ref KomaMRIBase.get_n_coils) to obtain its channel count:

x_axis = range(-FOV / 2, FOV / 2; length=recon_size[1])
y_axis = range(-FOV / 2, FOV / 2; length=recon_size[2])
x = vec([x for x in x_axis, _ in y_axis])
y = vec([y for _ in x_axis, y in y_axis])
z = zeros(length(x))
sensitivity_maps = reshape(
    get_sens(birdcage, x, y, z),
    recon_size..., 1, get_n_coils(birdcage),
);
sensitivity_maps = eltype(first(raw_coils_full.profiles).data).(sensitivity_maps); #hide

# [`BirdcageCoilSens`](@ref KomaMRIBase.BirdcageCoilSens) evaluates an analytic
# model. Measured maps can instead be supplied with
# [`ArbitraryCoilSens`](@ref KomaMRIBase.ArbitraryCoilSens); its
# [`get_sens`](@ref KomaMRIBase.get_sens) method interpolates those samples onto
# the same reconstruction grid.

# Four channels distributed around the array illustrate why it contains spatial information.
# The top row shows sensitivity magnitude and the bottom row shows the
# corresponding fully sampled coil image.

displayed_coils = 1:2:8 #hide
coil_image_scale = maximum(abs, coil_images_full[:, :, 1, 1, displayed_coils, 1]) #hide
sensitivity_plots = [ #hide
    image_panel( #hide
        abs.(sensitivity_maps[:, :, 1, coil]) ./ #hide
        maximum(abs, sensitivity_maps[:, :, 1, coil]), #hide
        "Sensitivity $coil"; zmax=1, #hide
    ) #hide
    for coil in displayed_coils #hide
] #hide
coil_image_plots = [ #hide
    image_panel( #hide
        abs.(coil_images_full[:, :, 1, 1, coil, 1]), #hide
        "Coil image $coil"; zmax=coil_image_scale, #hide
    ) #hide
    for coil in displayed_coils #hide
] #hide
p_coils = square_axes!( #hide
    [ #hide
        sensitivity_plots[1] sensitivity_plots[2] sensitivity_plots[3] sensitivity_plots[4] #hide
        coil_image_plots[1] coil_image_plots[2] coil_image_plots[3] coil_image_plots[4] #hide
    ]; #hide
    hide_ticks=true, height=380, #hide
) #hide
#md p_coils #hide
#jl display(p_coils)

# ## Recover the accelerated image with SENSE

# SENSE instead combines the unknown image with the complex sensitivity maps
# before sampling each coil's k-space trajectory. The inverse problem is solved
# jointly across all 16 channels:

sense_params = Dict{Symbol,Any}(
    :reco => "multiCoil",
    :reconSize => recon_size,
    :senseMaps => sensitivity_maps,
    :iterations => 20,
    :densityWeighting => false,
    :toeplitz => false,
)
sense_acc = reconstruction(acq_coils_acc, sense_params);

sense_acc_image = abs.(sense_acc[:, :, 1, 1, 1, 1]) #hide
p_sense_acc = image_panel( #hide
    sense_acc_image, "R = 3 SENSE"; zmax=maximum(sense_acc_image), #hide
) #hide
p_sos_acc_compare = image_panel( #hide
    sos_acc_image, "R = 3 SoS: aliased"; zmax=maximum(sos_acc_image), #hide
) #hide
p_sense = square_axes!( #hide
    [p_sos_acc_compare p_sense_acc]; #hide
    hide_ticks=true, height=280, #hide
) #hide
#md p_sense #hide
#jl display(p_sense)

# The ``R=3`` SoS image still folds because channel combination alone cannot
# separate overlapping pixels. SENSE uses the distinct complex channel
# sensitivities to remove that aliasing while preserving the shorter acquisition.

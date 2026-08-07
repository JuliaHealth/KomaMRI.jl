# Use Custom Coil Sensitivity Maps

This example estimates ESPIRiT sensitivity maps from a fully sampled reference
MRD, uses them in `ArbitraryCoilSens`, and compares direct and SENSE
reconstructions of measured and simulated `R=2` EPI data.

## Acquisition setup

The code below assumes `fully_sampled_mrd_file`, `accelerated_2x_mrd_file`,
`accelerated_2x_seq_file`, and `phantom_2D_3T_file` point to the four local
input files.

```julia
recon_size = (128, 128)
R = 2
seq = read_seq(accelerated_2x_seq_file)

number_of_shots = Int(seq.DEF["EpiShots"])
acquired_shot_indices =
    parse.(Int, split(seq.DEF["EpiAcquiredShots"], ',')) .- 1
acquired_line_indices = [
    line_index for shot_index in acquired_shot_indices for
    line_index in shot_index:number_of_shots:(recon_size[2] - 1)
]
@assert length(acquired_line_indices) == recon_size[2] ÷ R
```

`EpiShots` is the number of EPI interleaves. `EpiAcquiredShots` lists the
interleaves retained by the accelerated scan using one-based numbering, while
ISMRMRD phase-encoding indices start at zero. Subtracting one converts between
the two conventions.

For this sequence, there are 20 interleaves and the `R=2` scan keeps
`[0, 2, ..., 18]`. Interleaf `0` contains lines `[0, 20, ..., 120]`, interleaf
`2` contains `[2, 22, ..., 122]`, and so on. The comprehension preserves this
shot-by-shot acquisition order and produces 64 indices, which the assertion
checks against half of the 128 phase-encoding lines.

The raw files we have also begin with three navigator profiles that are not image
readouts. This helper discards them, keeps the first image average, and assigns
the corresponding phase-encoding index to every retained profile:

```julia
function select_profiles!(raw, line_indices)
    navigator_count = 3
    raw.profiles = raw.profiles[
        (navigator_count + 1):(navigator_count + length(line_indices))
    ]
    for (profile, line_index) in zip(raw.profiles, line_indices)
        profile.head.idx.kspace_encode_step_1 = UInt16(line_index)
    end
    nothing
end
```

## Sensitivity maps

```julia
raw_reference = RawAcquisitionData(ISMRMRDFile(fully_sampled_mrd_file))
raw_reference.profiles = raw_reference.profiles[4:(3 + recon_size[2])]

acq_reference = AcquisitionData(raw_reference)
acq_reference.traj[1].circular = false
sensitivity_maps = espirit(acq_reference, (6, 6), 30, recon_size; eigThresh_1=0.02, eigThresh_2=0.0,)
```

The reference MRD starts with three navigator profiles and then contains the
image profiles for each average. Selecting profiles `4:(3 + 128)` therefore
discards the navigators and keeps the first fully sampled image average.

`AcquisitionData` places those profiles on the reconstruction grid. The EPI
data cover a rectangular Cartesian grid, so `circular=false` prevents the
reconstruction from applying a circular shutter.

For this 12-coil reference, we use a conventional 6 × 6 kernel and a central
30 × 30 calibration region, which provides enough overlapping calibration patches
for a stable ESPIRiT estimate. `eigThresh_1=0.02` keeps calibration
singular vectors whose values are at least 2% of the largest one.
`eigThresh_2` controls the final spatial support mask.

## Measured acquisition

```julia
raw_measured = RawAcquisitionData(ISMRMRDFile(accelerated_2x_mrd_file))
select_profiles!(raw_measured, acquired_line_indices)
acq_measured = AcquisitionData(raw_measured)
acq_measured.traj[1].circular = false

direct_params = Dict{Symbol,Any}(:reco => "direct", :reconSize => recon_size,)
sense_params = Dict{Symbol,Any}(
    :reco => "multiCoil",
    :reconSize => recon_size,
    :senseMaps => sensitivity_maps,
    :iterations => 20,
    :densityWeighting => false,
    :toeplitz => false,
)
```

`select_profiles!` removes the three navigators and assigns the 64 acquired
phase-encoding lines to their correct positions in the `128 × 128` grid. The
direct reconstruction preserves the expected `R=2` aliasing in each coil,
whereas `multiCoil` uses the sensitivity maps to combine the 12 coils and
unfold the image.

```julia
magnitude_image(image) = abs.(Array(image[:, :, 1, 1, 1, 1]))
plot_reconstruction(image, title) = plot_image(image; title, zmin=0, zmax=quantile(vec(image), 0.995),)

measured_direct = reconstruction(acq_measured, direct_params)
measured_sense = reconstruction(acq_measured, sense_params)
measured_direct_image = reverse(magnitude_image(measured_direct); dims=(1, 2))
measured_sense_image = reverse(magnitude_image(measured_sense); dims=(1, 2))
p_measured_direct = plot_reconstruction(measured_direct_image, "Measured direct (R=2)")
p_measured_sense = plot_reconstruction(measured_sense_image, "Measured SENSE (R=2)")
```

`magnitude_image` selects the first slice, echo, coil, and repetition and takes
its magnitude. The direct image therefore shows the first receive coil, while
the SENSE result already contains one coil-combined image. The 99.5th percentile
sets a robust display maximum so isolated bright pixels do not reduce the
visible contrast. Both measured images are rotated by 180 degrees to match the
displayed phantom orientation.

Hover to inspect values, drag or scroll inside a figure to zoom, and double-click
to reset.

```@raw html
<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem;">
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/acquisition_direct.html" style="width:100%;height:520px;"></object>
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/acquisition_sense.html" style="width:100%;height:520px;"></object>
</div>
```

## Simulated acquisition

`ArbitraryCoilSens` expects maps in `(x, y, z, coil)` order and coordinates in
metres. The reference field of view is stored in millimetres, so it is converted
to metres before defining the map grid. The 2D maps are repeated at three
positions through the slice without changing the 12-coil dimension.

```julia
fov = Float32.(raw_reference.params["reconFOV"]) .* 1f-3
x = collect(LinRange(-fov[1] / 2, fov[1] / 2, recon_size[1]))
y = collect(LinRange(-fov[2] / 2, fov[2] / 2, recon_size[2]))
z = Float32[-fov[3] / 2, 0, fov[3] / 2]
receiver = ArbitraryCoilSens(x, y, z, repeat(sensitivity_maps, 1, 1, length(z), 1),)
```

`ArbitraryCoilSens` interpolates these maps at the phantom positions and returns
zero outside the supplied grid, so the grid must cover the complete phantom.
This is a triggered sequence, so the simulation uses `heart_rate=1`.

```julia
obj = read_phantom(phantom_2D_3T_file)
raw_simulated = simulate( obj, seq, Scanner(; receiver); physio=CardiacSignal(; heart_rate=1), verbose=false,)

simulated_mrd_file = "simulated_acquisition.mrd"
save(ISMRMRDFile(simulated_mrd_file), raw_simulated)
raw_simulated = RawAcquisitionData(ISMRMRDFile(simulated_mrd_file))
select_profiles!(raw_simulated, acquired_line_indices)
```

Saving and reopening the simulated data verifies the same ISMRMRD workflow used
for the measured scan. `select_profiles!` then removes the simulated navigators
and restores the accelerated phase-encoding indices.

```julia
acq_simulated = AcquisitionData(raw_simulated)
acq_simulated.traj[1].circular = false
simulated_direct = reconstruction(acq_simulated, direct_params)
simulated_sense = reconstruction(acq_simulated, sense_params)
simulated_direct_image = magnitude_image(simulated_direct)
simulated_sense_image = magnitude_image(simulated_sense)
p_simulated_direct = plot_reconstruction(simulated_direct_image, "Simulated MRD direct (R=2)")
p_simulated_sense = plot_reconstruction(simulated_sense_image, "Simulated SENSE (R=2)")
```

The same direct and SENSE settings used for the measured acquisition are reused
here, making the two pairs of results directly comparable. The simulated images
already match the phantom orientation and therefore are not rotated.

```@raw html
<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem;">
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/simulated_mrd_direct.html" style="width:100%;height:520px;"></object>
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/simulated_mrd_sense.html" style="width:100%;height:520px;"></object>
</div>
```

## Cardiac MOLLI

The fully sampled central calibration region from `SA_test` provides 30
RSS-normalized coil sensitivity estimates. Use the slider to select a channel.

```@raw html
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/sa_test_coil_sensitivities.html" style="width:100%;height:420px;"></object>
```

The `SA_test` maps are used directly for the accelerated `R=2` MOLLI scan.
Direct reconstruction retains the two-fold aliasing; `multiCoil` unfolds it.

```@raw html
<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem;">
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/cardiac_molli_direct.html" style="width:100%;height:520px;"></object>
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/cardiac_molli_multi_coil.html" style="width:100%;height:520px;"></object>
</div>
```

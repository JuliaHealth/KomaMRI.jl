# KomaMRIFiles

```@meta
CurrentModule = KomaMRIFiles
```

## Sequence

```@docs
read_seq
read_seq_data
write_seq
write_seq_data
PulseqSequenceData
```

## Phantom

```@docs
read_phantom_jemris
read_phantom_MRiLab
read_phantom
write_phantom
```

## Scanner

Create or load an HDF5 `.sys` file:

```julia
write_scanner(sys, "scanner.sys")
sys = read_scanner("scanner.sys")
```

```@docs
read_scanner
write_scanner
```

## MATLAB exports

The following formats are used by KomaMRI's UI and CLI exports.
UI individual exports use `data_raw.mat`, `data_sequence.mat`, and `data_image.mat`;
`All` uses `raw.mat`, `seq_sequence.mat`, and `image.mat`.
CLI output filenames are user-specified.

**Format change:** `raw` is now a struct, replacing the old numeric `[time, signal]`
matrix. Label-aware image exports now contain a `reconstruction` struct instead of the
previous separate `image`, `labels`, `source_labels`, and `recon_policy` variables.
Plain-array image exports still contain the original `image` variable.

### Raw data

`raw.params` preserves the MRD metadata. `raw.profiles` is a cell array, with one struct
per acquired readout:

| Field | Contents |
| --- | --- |
| `head` | Full MRD acquisition header, including flags, timestamps, discard counts, and `idx` label counters |
| `traj` | Stored trajectory, coordinates × samples; normalization is unchanged |
| `data` | Complex signal, samples × receive coils |
| `role` | Acquisition role derived from the MRD flags, such as `imaging` or `navigator` |

Profiles retain their own sample counts and trajectory dimensions, so 1D navigators,
2D images, and 3D acquisitions can coexist. No coils are discarded, no `Inf` separators
are inserted, and timestamps are not synthesized. Koma's trajectory normalization scale
is retained in `raw.params.userParameters.KomaTrajectoryScale` when available.

Raw data also exports when `userParameters` is absent. When present, those parameters
are additionally saved in `sim_params.mat`. Dictionary keys replace `Δ` with `d` for
MATLAB compatibility, for example `Δt_rf` becomes `dt_rf`.

```matlab
s = load('data_raw.mat');
p = s.raw.profiles{1};
signal = p.data(:, 1);       % All samples from receive coil 1
k = p.traj;                 % Original stored coordinates for this readout
slice_label = p.head.idx.slice;
```

### Sequence data

`sequence` retains the waveform fields `Gx`, `Gy`, `Gz`, `RF_AM`, `RF_FM`, and `ADCS`,
and adds `definitions` and `adc`. Each `sequence.adc` vector has one entry per ADC event:
`block` is the one-based Julia sequence-block index, `num_samples` is its sample count,
and `labels` contains named vectors of all accumulated ADC label values. Label values
are preserved, including zero; they are not converted to MATLAB indices.

```matlab
s = load('data_sequence.mat');
adc_blocks = s.sequence.adc.block;
repetitions = s.sequence.adc.labels.REP;
```

### Reconstructed images

`reconstruction.policy` records the reconstruction policy. Each cell in
`reconstruction.images` contains one reconstruction batch with `data`, `size`, `labels`,
and `source`. `labels` identifies the batch; `source` retains its original
acquisition counters. The six data dimensions are **x, y, z, echo, coil, repetition**.
The explicit six-element `size` vector preserves trailing singleton dimensions that
MATLAB does not display.

```matlab
s = load('data_image.mat');
batch = s.reconstruction.images{1};
I = reshape(batch.data, batch.size(:).');
echo = 1; coil = 1; repetition = 1; z = 1;
plane = I(:, :, z, echo, coil, repetition);
if batch.size(2) == 1 && batch.size(3) == 1
    plot(abs(plane(:, 1)));       % 1D reconstruction
else
    imagesc(abs(plane)); axis image; % 2D image or one plane of a 3D volume
end
volume = I(:, :, :, echo, coil, repetition);
```

Select a different batch for separate acquisition labels, such as `SLC` or `REP`.
Select `z` within a batch for a reconstructed 3D partition. MATLAB array indices are
one-based, independently of the stored label values. Explicit indexing above avoids
collapsing coil or echo dimensions with `squeeze`.

Profile and image collections remain cell arrays even with one element, for a
consistent format across supported MAT.jl versions 0.10, 0.11, and 0.12.

### Scanner data

`scanner` contains `limits`, `gradient`, `receiver`, and `transmitter` structs,
each with its model `type` and fields. For example, `scanner.limits.B0` is the
main field strength, and `scanner.receiver.coil_sens` stores complex maps for
`ArbitraryCoilSens`. This replaces the previous flat limits-only struct.

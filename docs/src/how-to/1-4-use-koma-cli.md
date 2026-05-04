# Use Koma from the Command Line

KomaMRI provides the `koma` Julia app for opening the UI or running simulations from a terminal.

## Installing the App

Install KomaMRI with `pkg> add KomaMRI`, then install the `koma` app (Julia 1.12+):
```julia-repl
pkg> app add KomaMRI # Julia 1.12+
pkg> # app dev /path/to/KomaMRI.jl # Local checkout, Julia 1.12+
```

The executable is installed in `~/.julia/bin`. If `koma` is not found, add that directory to your `PATH`.
The app command is named `koma`; Julia app names are executable names, not registered package names.

## Basic Workflow & Input Defaults

With no arguments, `koma` opens the graphical interface with default inputs:
```bash
koma
```

Add input files without output paths to preload the UI before it opens:
```bash
koma -i epi.seq brain.phantom
```

Input files can be passed in any order. KomaMRI identifies them by extension: `.seq`, `.phantom` or `.h5`.
A scanner `.sys` file is also accepted, but it does nothing for now.

## Running Simulation

Add a simulation output path with `-o` or `--outputs` to run without opening the UI:
```bash
koma -i epi.seq brain.phantom -o raw.mrd
```

Simulation output can be `.mrd` or `.mat`, selected from the output filename extension.

## Reconstructing Image using MRIReco

Pass a second output path to reconstruct after simulation. Reconstruction output is currently `.mat`:
```bash
koma -i epi.seq brain.phantom -o raw.mrd image.mat
```

Use `_` to skip saving the simulation output:
```bash
koma -i epi.seq brain.phantom -o _ image.mat
```

## Simulation and Reconstruction Parameters

Simulation and reconstruction parameters can be passed as repeated `KEY=VALUE` pairs:
```bash
koma -i epi.seq brain.phantom -s sim_method=BlochMagnus4 -r reco=direct -o raw.mrd
```

The long forms are `--sim-param` and `--recon-param`.

## Backend and Threading

Simulations run on multi-threaded CPU by default. GPU backends can be loaded by name:
```bash
koma -b Metal -i epi.seq brain.phantom -o raw.mrd
```

The selected backend package must be installed in your default Julia environment.

The app starts with CPU threading enabled automatically. To choose a thread count, pass Julia flags before `--` and KomaMRI arguments after it:
```bash
koma --threads=8 -- -i epi.seq brain.phantom -o raw.mrd
```

## Default App Configuration

`koma` can read CLI defaults from Julia's `LocalPreferences.toml` in the app environment. Put the file in the path that matches how the app was installed:

- Installed with `pkg> app add KomaMRI` (Julia 1.12+): `~/.julia/environments/apps/KomaMRI/LocalPreferences.toml`
- Installed with `pkg> app dev /path/to/KomaMRI.jl` (Julia 1.12+): `/path/to/KomaMRI.jl/LocalPreferences.toml`

```toml
[KomaMRI.koma]
backend = "CPU" # CPU, CUDA, Metal, AMDGPU, or oneAPI

[KomaMRI.koma.inputs]
sequence = "mysequence.seq"
phantom = "myphantom.phantom"

[KomaMRI.koma.outputs]
# rawdata = "raw.mrd"
# image = "image.mat"

[KomaMRI.koma.sim_params]
# Numerical accuracy
sim_method = "BlochMagnus4"
precision = "f32"
"Δt" = 1e-3
"Δt_rf" = 5e-5
# Memory
max_block_length = 512
max_rf_block_length = inf
# CPU multi-threaded
Nthreads = 8
# GPU kernel launch config
gpu_groupsize_precession = 256
gpu_groupsize_excitation = 256
```

CLI arguments override these defaults.

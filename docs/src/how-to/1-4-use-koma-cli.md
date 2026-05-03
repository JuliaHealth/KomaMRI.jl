# Use Koma from the Command Line

KomaMRI provides the `koma` Julia app for opening the UI or running simulations from a terminal.

## Installing the App

Install KomaMRI with `pkg> add KomaMRI`, then install the `koma` app (Julia 1.12+):
```julia-repl
julia> import Pkg
julia> Pkg.Apps.add("KomaMRI") # Installs koma app
julia> # Pkg.Apps.develop(path="/path/to/KomaMRI.jl") # Local checkout
```

The executable is installed in `~/.julia/bin`. If `koma` is not found, add that directory to your `PATH`.
The app command is named `koma`; Julia app names are executable names, not registered package names.

## Basic Workflow

With no arguments, `koma` opens the graphical interface with default inputs:
```bash
koma
```

Add input files without output paths to preload the UI before it opens:
```bash
koma -i seq.seq obj.phantom sys.sys
```

Input files can be passed in any order. KomaMRI identifies them by extension: `.seq`, `.phantom` or `.h5`, and `.sys`.

## Loading Simulation Inputs

Use `-i` or `--inputs` to pass a sequence, phantom, and scanner file. Missing inputs use the same defaults as the UI:
```bash
koma -i obj.phantom seq.seq sys.sys
```

Scanner `.sys` files are accepted but ignored for now.

## Running Simulation

Add a simulation output path with `-o` or `--outputs` to run without opening the UI:
```bash
koma -i obj.phantom seq.seq sys.sys -o raw.mrd
```

Simulation output can be `.mrd` or `.mat`, selected from the output filename extension.

## Reconstructing Image using MRIReco

Pass a second output path to reconstruct after simulation. Reconstruction output is currently `.mat`:
```bash
koma -i seq.seq obj.phantom sys.sys -o raw.mrd image.mat
```

Use `_` to skip saving the simulation output:
```bash
koma -i seq.seq obj.phantom sys.sys -o _ image.mat
```

## Simulation and Reconstruction Parameters

Simulation and reconstruction parameters can be passed as repeated `KEY=VALUE` pairs:
```bash
koma -s gpu=false -r reco=direct -o raw.mrd
```

The long forms are `--sim-param` and `--recon-param`.

## Backend and Threading

Simulations run on multi-threaded CPU by default. GPU backends can be loaded by name:
```bash
koma -b Metal -i seq.seq obj.phantom sys.sys -o raw.mrd
```

The selected backend package must be installed in your default Julia environment.

The app starts with CPU threading enabled automatically. To choose a thread count, pass Julia flags before `--` and KomaMRI arguments after it:
```bash
koma --threads=8 -- -i seq.seq obj.phantom sys.sys -o raw.mrd
```

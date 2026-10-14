# Use Koma in Julia Scripts

Run the same brain/EPI example as the [UI](1-1-use-koma-ui.md#Basic-Workflow)
and [notebook](1-2-use-koma-notebooks.md#Basic-Workflow) guides from a Julia script or REPL.

## Basic Workflow
(You can also go to [analog steps using UI](1-1-use-koma-ui.md#Basic-Workflow))

As a general overview, remember the following workflow steps when using KomaMRI:

* Loading Simulation Inputs: **Scanner**, **Phantom**, **Sequence**
* Running Simulation
* Reconstructing Image using **MRIReco**

Let's replicate these previous steps in a **Julia** script. You will obtain the following code, which you can copy and paste into the **Julia REPL**:
Install `MAT` with `pkg> add MAT` if needed for the export step.

```julia
# Import the package
using KomaMRI, MAT

# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
obj.Δw .= 0 # Match the UI's default phantom
seq = PulseDesigner.EPI_example(; sys)

# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; sim_params)

# Perform reconstruction to get the image
rec_params = Dict{Symbol,Any}(:reco => "direct")
result = reconstruct_with_labels(raw; rec_params)
Nx, Ny = raw.params["reconSize"][1:2]
image = reshape(first(result.images).image, Nx, Ny)

# Export raw data and the image
save(ISMRMRDFile(joinpath(tempdir(), "data_raw.mrd")), raw)
MAT.matwrite(joinpath(tempdir(), "data_image.mat"), Dict("image" => image))
```

Let's go through this script step by step.

## Loading Simulation Inputs
(You can also go to [analog steps using UI](1-1-use-koma-ui.md#Loading-Simulation-Inputs))

The inputs of the simulation are created in the following part of the script: 
```julia
# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
obj.Δw .= 0
seq = PulseDesigner.EPI_example(; sys)
```

```@setup script_inputs
using KomaMRI
sys = Scanner()
obj = brain_phantom2D()
obj.Δw .= 0
seq = PulseDesigner.EPI_example(; sys)
sim_params = KomaMRICore.default_sim_params()
```

### Scanner

[`Scanner`](@ref) groups `limits`, `gradient`, `receiver`, and `transmitter`.
Inspect the hardware limits with `sys.limits`:

```@repl script_inputs
sys.limits
```

### Phantom

The Phantom struct created in this example represents a slice of a brain. To create it, we use the function `brain_phantom2D`, which is part of **KomaMRIBase**. While **KomaMRI** provides some phantom examples for experimentation, you may also want to create your custom **Phantom** struct tailored to your specific requirements.

The **Phantom** struct contains MRI parameters related to the magnetization properties of an object. These parameters include magnetization positions, proton density, relaxation times, off-resonance, among others.
For more information about Koma's Phantom and what it can do, as well as how to store and share it, check out the [Phantom](../explanation/gen-1-phantom.md) section.
Inspect the first spin's parameters in the **Julia REPL**:
```@repl script_inputs
obj[1]
```

You can also visualize the **Phantom** struct using the [`plot_phantom_map`](@ref) function:
```julia-repl
julia> plot_phantom_map(obj, :ρ)
```
```@raw html
<object type="text/html" data="../assets/phantom-rho.html" style="width:100%; height:620px;"></object>
```

To utilize test phantoms included with **KomaMRI**, navigate to the "examples" folder and use the [`read_phantom_jemris`](@ref)  function to read a phantom in `.h5` format. The following steps outline how to do this in **Julia**:
```julia-repl
julia> path_koma = dirname(dirname(pathof(KomaMRI)))
julia> path_sphere = joinpath(path_koma, "examples", "2.phantoms", "sphere_chemical_shift.h5")
julia> sphere = read_phantom_jemris(path_sphere)
julia> plot_phantom_map(sphere, :T2)
```
```@raw html
<object type="text/html" data="../assets/phantom-T2-circle.html" style="width:100%; height:620px;"></object>
```

This explores `sphere` without replacing the brain phantom `obj` used below.

### Sequence

The **Sequence** struct in the example represents one of the most basic MRI sequences. It excites the object with a 90° RF pulse and then uses EPI gradients to fill the k-space in a "square" manner. While you may want to create your sequences for experiments, you can always use some of the examples already available in **KomaMRI**.

In MRI, the sequence must be carefully designed with precise timing to obtain an image. It includes subcomponents such as gradients, radio-frequency excitation signals, and sample acquisition. For more information on constructing a **Sequence** struct, refer to the [Sequence](../explanation/4-sequence.md) section.

You can view general information about a **Sequence** struct by displaying it in the **Julia REPL**:
```@repl script_inputs
seq
```

For more precise timing checks, you can use the [`plot_seq`](@ref) function:
```julia-repl
julia> plot_seq(seq; range=[0 30])
```
```@raw html
<object type="text/html" data="../assets/plot-seq-epi.html" style="width:100%; height:420px;"></object>
```

It is important to consider how the sequence traverses through k-space. The [`plot_kspace`](@ref) function does precisely that:
```julia-repl
julia> plot_kspace(seq)
```
```@raw html
<object type="text/html" data="../assets/kspace-epi.html" style="width:100%; height:420px;"></object>
```

For new sequences, use **PulseDesigner** constructors such as
[`PulseDesigner.make_block_pulse`](@ref), [`PulseDesigner.make_sinc_pulse`](@ref),
[`PulseDesigner.make_trapezoid`](@ref), and [`PulseDesigner.make_adc`](@ref).
See [Build Sequences with `@addblock`](3-create-your-own-sequence.md).

**KomaMRI** is also compatible with [Pulseq](https://pulseq.github.io/). The package installation includes some `.seq` files in **Pulseq** format, which can be read and used as a **Sequence** struct. Here's how to read a spiral **Pulseq** file stored in the "examples" folder of **KomaMRI**:
```julia-repl
julia> path_koma = dirname(dirname(pathof(KomaMRI)))
julia> path_spiral = joinpath(path_koma, "examples", "1.sequences", "spiral.seq")
julia> spiral = read_seq(path_spiral)
julia> plot_seq(spiral)
julia> plot_kspace(spiral)
```
```@raw html
<object type="text/html" data="../assets/seq-spiral-pulseq-time.html" style="width:50%; height:420px;"></object><object type="text/html" data="../assets/seq-spiral-pulseq-kspace.html" style="width:50%; height:420px;"></object>
```

This explores `spiral` without replacing the EPI sequence `seq` used below.

## Running Simulation
(You can also go to [analog steps using UI](1-1-use-koma-ui.md#Running-Simulation))

The following lines in the example script configure and perform the simulation:
```julia
# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; sim_params)
```

### Simulation Parameters

To perform simulations, **KomaMRI** requires certain parameters. You can use the default parameters for testing, but you also have the option to customize specific simulation aspects. In the example, we use the [`KomaMRICore.default_sim_params`](@ref) function to create a dictionary with default simulation parameters. You can view the keys that can be modified by displaying the `sim_params` variable:
```@repl script_inputs
sim_params
```

`"sampling_rule"` controls sequence sampling. The compatibility keys `"Δt"` and
`"Δt_rf"` specify maximum time steps, not scanner raster times. `"return_type"`
defaults to MRD raw data; `"mat"` returns a signal array. `"gpu"` selects GPU use,
and `"precision"` sets floating-point precision. See
[Simulation Parameters](../explanation/7-simulation.md).


### Raw Signal

The simulation is performed using the [`simulate`](@ref) function, which requires three arguments: a **Scanner** struct, a **Phantom** struct, and a **Sequence** struct. Optionally, you can include the keyword argument `sim_params` if you wish to use custom simulation parameters.

In the example, we can see that the output of the simulation is a special struct:
```@repl script_inputs
raw = simulate(obj, seq, sys; sim_params=merge(sim_params, Dict("gpu" => false)), verbose=false); # hide
typeof(raw)
raw
```

You can plot the simulation result with the [`plot_signal`](@ref) function like so:
```julia-repl
julia> plot_signal(raw)
```
```@raw html
<object type="text/html" data="../assets/raw-epi-brain-default.html" style="width:100%; height:420px;"></object>
```


## Reconstructing Image using MRIReco
(You can also go to [analog steps using UI](1-1-use-koma-ui.md#Reconstructing-Image-using-MRIReco))

Use [`reconstruct_with_labels`](@ref), the same MRIReco-based reconstruction
used by the UI. This example contains one 2D image:
```julia
rec_params = Dict{Symbol,Any}(:reco => "direct")
result = reconstruct_with_labels(raw; rec_params)
Nx, Ny = raw.params["reconSize"][1:2]
image = reshape(first(result.images).image, Nx, Ny)
```

See the [reconstruction reference](../reference/6-koma-mri.md#Reconstruction)
for acquisition labels, coil combination, and reconstruction policies.

Display the image with the same magnitude scaling as the UI:
```julia-repl
julia> plot_image(image * (Nx * Ny))
```
```@raw html
<div style="text-align:center;"><object type="text/html" data="../assets/image-default-brain.html" style="width:100%; height:620px;"></object></div>
```


## Exporting Results

Same outputs: [UI](1-1-use-koma-ui.md#Exporting-results).

Save the existing raw data as MRD and the reconstructed image as MAT:
```julia
using MAT
save(ISMRMRDFile(joinpath(tempdir(), "data_raw.mrd")), raw)
MAT.matwrite(joinpath(tempdir(), "data_image.mat"), Dict("image" => image))
```

The files are in `tempdir()`; replace that path with your output directory.

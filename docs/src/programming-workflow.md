# Julia Scripts

You should already be familiar with the  [Graphical User Interface](ui-details.md) of **KomaMRI**. However, you can also use this package directly from the **Julia REPL** or write your own Julia scripts. This allows you to unlock the full potential of KomaMRI, enabling you to utilize more of its functionalities and even test your own MRI ideas.

This section demonstrates a basic workflow with **KomaMRI** through writing your own scripts or entering commands directly into the **Julia REPL**. Let's begin.

## Basic Workflow
(You can also go to [analog steps using UI](ui-details.md#Basic-Workflow))

As a general overview, remember the following workflow steps when using KomaMRI:

* Loading Simulation Inputs: **Scanner**, **Phantom**, **Sequence**
* Running Simulation
* Reconstructing Image using **MRIReco**

Let's replicate these previous steps in a **Julia** script. You will obtain the following code, which you can copy and paste into the **Julia REPL**:
```julia
# Import the package
using KomaMRI

# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
seq = PulseDesigner.EPI_example()

# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; sim_params)

# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
    Nx, Ny = raw.params["reconSize"][1:2]
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    image3d  = reshape(rec.data, Nx, Ny, :)
    image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
    return image2d
end

# Perform reconstruction to get the image
image = reconstruct_2d_image(raw)
```

Let's go through this script step by step.

## Loading Simulation Inputs
(You can also go to [analog steps using UI](ui-details.md#Loading-Simulation-Inputs))

The inputs of the simulation are created in the following part of the script: 
```julia
# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
seq = PulseDesigner.EPI_example()
```

### Scanner

The previously created **Scanner** struct contains default parameters. In your initial simulations, you will likely use this default struct without making any modifications. You can view all the parameters by displaying the struct variable in the **Julia REPL**. The Scanner's parameters include hardware limitations such as the main magnetic field, maximum gradient values, minimum raster times, and more. You may want to adjust these values for your future custom simulations.
```julia-repl
julia> sys
Scanner
  B0: Float64 1.5
  B1: Float64 1.0e-5
  Gmax: Float64 0.06
  Smax: Int64 500
  ADC_Δt: Float64 2.0e-6
  seq_Δt: Float64 1.0e-5
  GR_Δt: Float64 1.0e-5
  RF_Δt: Float64 1.0e-6
  RF_ring_down_T: Float64 2.0e-5
  RF_dead_time_T: Float64 0.0001
  ADC_dead_time_T: Float64 1.0e-5
```

### Phantom

The Phantom struct created in this example represents a slice of a brain. To create it, we use the function **brain\_phantom2D()**, which is part of the subdependency **KomaMRICore**. While **KomaMRI** provides some phantom examples for experimentation, you may also want to create your custom **Phantom** struct tailored to your specific requirements.

The **Phantom** struct contains MRI parameters related to the magnetization properties of an object. These parameters include magnetization positions, proton density, relaxation times, off-resonance, among others. To view all the keys and values of the object, you can do so in the **Julia REPL** as follows:
```julia-repl
julia> obj
Phantom{Float64}
  name: String "brain2D_axial"
  x: Array{Float64}((6506,)) [-0.084, -0.084,  …  0.086, 0.086]
  y: Array{Float64}((6506,)) [-0.03, -0.028,  …  0.0, 0.002]
  z: Array{Float64}((6506,)) [-0.0, -0.0,  …  0.0, 0.0]
  ρ: Array{Float64}((6506,)) [0.7, 0.7,  …  0.7, 0.7]
  T1: Array{Float64}((6506,)) [0.569, 0.569,  …  0.569, 0.569]
  T2: Array{Float64}((6506,)) [0.329, 0.329,  …  0.329, 0.329]
  T2s: Array{Float64}((6506,)) [0.058, 0.058,  …  0.058, 0.058]
  Δw: Array{Float64}((6506,)) [-0.0, -0.0,  …  -0.0, -0.0]
  Dλ1: Array{Float64}((6506,)) [0.0, 0.0,  …  0.0, 0.0]
  Dλ2: Array{Float64}((6506,)) [0.0, 0.0,  …  0.0, 0.0]
  Dθ: Array{Float64}((6506,)) [0.0, 0.0,  …  0.0, 0.0]
...
```
As you can see, attributes of the **Phantom** struct are vectors representing object properties, with each element holding a value associated with a single magnetization.

You can also visualize the **Phantom** struct using the **plot\_phantom\_map()** function, which is part of the **KomaMRIPlots** subdependency. This function plots the magnitude of a property for each magnetization at a specific spatial position. You can observe properties such as proton density and relaxation times, so feel free to replace the **:ρ** symbol with another property of the phantom in the example below:
```julia-repl
julia> plot_phantom_map(obj, :ρ)
```
```@raw html
<object type="text/html" data="../assets/phantom-rho.html" style="width:100%; height:620px;"></object>
```

To utilize test phantoms included with **KomaMRI**, navigate to the "examples" folder and use the **read\_phantom\_jemris()** function to read a phantom in **.h5** format. The following steps outline how to do this in **Julia**:
```julia-repl
julia> path_koma = dirname(dirname(pathof(KomaMRI)))
julia> path_sphere = joinpath(path_koma, "examples", "2.phantoms", "sphere_chemical_shift.h5")
julia> sphere = read_phantom_jemris(path_sphere)
julia> plot_phantom_map(sphere, :T2)
```
```@raw html
<object type="text/html" data="../assets/phantom-T2-circle.html" style="width:100%; height:620px;"></object>
```

### Sequence

The **Sequence** struct in the example represents one of the most basic MRI sequences. It excites the object with a 90° RF pulse and then uses EPI gradients to fill the k-space in a "square" manner. While you may want to create your sequences for experiments, you can always use some of the examples already available in **KomaMRI**.

In MRI, the sequence must be carefully designed with precise timing to obtain an image. It includes subcomponents such as gradients, radio-frequency excitation signals, and sample acquisition. For more information on constructing a **Sequence** struct, refer to the [Sequence](sequence.md) section.

You can view general information about a **Sequence** struct by displaying it in the **Julia REPL**:
```julia-repl
julia> seq
Sequence[ τ = 62.846 ms | blocks: 204 | ADC: 101 | GR: 205 | RF: 1 | DEF: 5 ]
```

For more precise timing checks, you can use the **plot\_seq()** function:
```julia-repl
julia> plot_seq(seq; range=[0 30])
```
```@raw html
<object type="text/html" data="../assets/plot-seq-epi.html" style="width:100%; height:420px;"></object>
```

It is important to consider how the sequence traverses through k-space. The **get\_kspace()** function does precisely that:
```julia-repl
julia> plot_kspace(seq)
```
```@raw html
<object type="text/html" data="../assets/kspace-epi.html" style="width:100%; height:420px;"></object>
```

Additionally, there are helpful sequence construction functions within a submodule of **KomaMRI** called **PulseDesigner**. These functions include **RF\_hard()**, **RF\_sinc()**, **EPI()**, **radial\_base()** and **spiral\_base()**. For more details on how to use them, refer to the [API documentation](api.md).

**KomaMRI** is also compatible with [Pulseq](https://pulseq.github.io/). The package installation includes some **.seq** files in **Pulseq** format, which can be read and used as a **Sequence** struct. Here's how to read a spiral **Pulseq** file stored in the "examples" folder of **KomaMRI**:
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

## Running Simulation
(You can also go to [analog steps using UI](ui-details.md#Running-Simulation))

The following lines in the example script configure and perform the simulation:
```julia
# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; sim_params)
```

### Simulation Parameters

To perform simulations, **KomaMRI** requires certain parameters. You can use the default parameters for testing, but you also have the option to customize specific simulation aspects. In the example, we use the **default\_sim\_params()** function to create a dictionary with default simulation parameters. You can view the keys that can be modified by displaying the `sim_params` variable:
```julia-repl
julia> sim_params
Dict{String, Any} with 9 entries:
  "return_type" => "raw"
  "Nblocks"     => 20
  "gpu"         => true
  "Nthreads"    => 1
  "gpu_device"  => 0
  "sim_method"  => Bloch()
  "precision"   => "f32"
  "Δt"          => 0.001
  "Δt_rf"       => 5.0e-5
```

All of these parameters merit special attention. Here is a description of each:

* `"return_type"`: defines the output of the **simulation()** function. Possible values are `"raw"`, `"mat"`, and `"state"`, corresponding to outputting a **MRIReco** `RawAcquisitionData`, the signal values, and the last magnetization state of the simulation, respectively.
* `"sim_method"`: defines the type of simulation. The default value is **Bloch()**, but you can alternatively use the **BlochDict()** simulation method. Moreover, you have the flexibility to create your own methods without altering the **KomaMRI** source code; for further details, refer to the [Simulation Method Extensibility section](mri-theory.md#Simulation-Method-Extensibility).
* `"Δt"`: raster time for gradients.
* `"Δt_rf"`: raster time for RFs.
* `"precision"`: defines the floating-point simulation precision. You can choose between `"f32"` and `"f64"` to use `Float32` and `Float64` primitive types, respectively. It's important to note that, especially for GPU operations, using `"f32"` is generally much faster.
* `"Nblocks"` divides the simulation into a specified number of time blocks. This parameter is designed to conserve RAM resources, as **KomaMRI** computes a series of simulations consecutively, each with the specified number of blocks determined by the value of `"Nblocks"`.
* `"Nthreads"`: divides the **Phantom** into a specified number of threads. Because spins are modeled independently of each other, **KomaMRI** can solve simulations in parallel threads, speeding up the execution time.
* `"gpu"`: is a boolean that determines whether to use GPU or CPU hardware resources, as long as they are available on the host computer.
* `"gpu_device"`: sets the index ID of the available GPU in the host computer.

For instance, if you want to perform a simulation on the CPU with float64 precision using the **BlochDict()** method, you can do so like this:
```julia
# Set non-default simulation parameters and run simulation
sim_params = KomaMRICore.default_sim_params() 
sim_params["gpu"] = false
sim_params["precision"] = "f64"
sim_params["sim_method"] = BlochDict()
raw = simulate(obj, seq, sys; sim_params)
```


### Raw Signal

The simulation is performed using the **simulate()** function, which requires three arguments: a **Scanner** struct, a **Phantom** struct, and a **Sequence** struct. Optionally, you can include the keyword argument `sim_params` if you wish to use custom simulation parameters.

In the example, we can see that the output of the simulation is a special struct:
```julia-repl
julia> typeof(raw)
RawAcquisitionData

julia> raw
RawAcquisitionData[SeqName: epi | 101 Profile(s) of 101×1]
```

You can plot the simulation result with the **plot\_signal()** function like so:
```julia-repl
julia> plot_signal(raw)
```
```@raw html
<object type="text/html" data="../assets/raw-epi-brain-default.html" style="width:100%; height:420px;"></object>
```


## Reconstructing Image using MRIReco
(You can also go to [analog steps using UI](ui-details.md#Reconstructing-Image-using-MRIReco))

**KomaMRI** does not handle reconstruction; instead, you should utilize the **MRIReco** package to generate an image. For convenience, when you install **KomaMRI**, you also install **MRIReco**, allowing you to access functions from that package. You should pay special attention to the **RawAcquisitionData** and **AcquisitionData** structs, as well as the **reconstruction()** function.

In the example below, we define an auxiliary function, **reconstruct\_2d\_image()**, which takes a raw signal struct, **RawAcquisitionData**, as input and returns a 2D Array representing an image. Within this function, we create an **AcquisitionData** struct and set some reconstruction parameters, which serve as inputs for the **reconstruction()** function. The latter function is responsible for the image generation process.
```julia
# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
    Nx, Ny = raw.params["reconSize"][1:2]
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    image3d  = reshape(rec.data, Nx, Ny, :)
    image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
    return image2d
end

# Perform reconstruction to get the image
image = reconstruct_2d_image(raw)
```

If you need more information about how to use the **AcquisitionData** and the how to fill the reconstruction parameters, you need to visit the [MRIReco webpage](https://github.com/MagneticResonanceImaging/MRIReco.jl)).

To display the image, you can use the **plot_\image()** function which is part of the **KomaMRIPlots** subpackage:
```julia-repl
julia> plot_image(image)
```
```@raw html
<center><object type="text/html" data="../assets/image-default-brain.html" style="width:100%; height:620px;"></object></center>
```


## Exporting Results to .mat File
(You can also go to [analog steps using UI](ui-details.md#Exporting-Results-to-.mat-File))

Many people in the MRI community uses MATLAB, probably you are one of them and you want to process the raw signal in the MATLAB environment after simulation is done with **KomaMRI**. Here we show you an example of how to save a **.mat** file with the information of the raw signal thank to the help of the **MAT** package:

Many people in the MRI community use **MATLAB**; you might be one of them and may want to process the **Raw Signal** in the **MATLAB** environment after simulation is completed with **KomaMRI**. Here, we provide an example of how to save a **.mat** file containing the  **Raw Signal** information using the **MAT** package.
```julia
# Use the MAT package
using MAT

# Perform simulation to return an Array type
sim_params["return_type"] = "mat"
raw = simulate(obj, seq, sys; sim_params)

# Save the .mat file in the temp directory
matwrite(joinpath(tempdir(), "koma-raw.mat"), Dict("raw" => raw))
```

Note that we need to simulate to return an array type (not the default **RawAcquisitionData**), and then we utilize the **matwrite()** function to save a file named "koma-raw.mat" in your computer's temporary directory. Now, you can navigate to your temporary directory (which you can find by displaying the result of **tempdir()** in the **Julia REPL**) and locate the "koma-raw.mat" file.

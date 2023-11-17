# Ways of using Koma

**KomaMRI** can be used in different environments, depending on the degree of flexibility you desire. If you lack prior programming knowledge, we recommend starting with the [User Interface](#User-Interface). If you seek the full range of flexibility that **KomaMRI** offers, programming with [Scripts](#Scripts) is likely your preference. Alternatively, you can utilize the programming environment provided by [Notebooks](#Notebooks). Detailed explanations for each method of using **KomaMRI** will be provided in the following sections.

## User Interface

This section explains how to use the user interface of the **KomaMRI** package and the internal processes during interaction.

### Basic Workflow
(You can also go to [analog steps using Scripts](programming-workflow.md#Basic-Workflow))

As a general overview, remember the following workflow steps when using KomaMRI:

* Loading Simulation Inputs: **Scanner**, **Phantom**, **Sequence**
* Running Simulation
* Reconstructing Image using **MRIReco**

In the following subsections, we will cover all the mentioned steps. First, open the **Julia REPL** and enter the following commands to include the **KomaMRI** package and launch the user interface:
```julia-repl
julia> using KomaMRI

julia> KomaUI()
```
```@raw html
<p align="center"><img width="90%" src="../assets/gui-dashboard.png"/></p>
```

### Loading Simulation Inputs
(You can also go to [analog steps using Scripts](programming-workflow.md#Loading-Simulation-Inputs))

The user interface has preloaded certain inputs into RAM, including the **Scanner**, **Phantom**, and **Sequence** structs. In the following subsections, we will demonstrate how to visualize these inputs.

#### Scanner

You can visualize the preloaded **Scanner** struct by clicking on the `Scanner` dropdown and then pressing the `View Scanner` button. The **Scanner** struct contains hardware-related information, such as the main magnetic field's magnitude:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-scanner-view.png"/></p>
```

#### Phantom

To see the phantom already stored in RAM, simply click on the `Phantom` dropdown an then press the `View Phantom` button. The preloaded phantom is a slice of a brain:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-phantom-view.png"/></p>
```

It is also possible to load `.h5` phantom files. The **KomaMRI.jl** has some examples stored at `~/.julia/packages/KomaMRI/<id-string>/examples/2.phantoms/`. For instance, let's load the `sphere_chemical_shift.h5` file:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-phantom-load.png"/></p>
```

Note that you can select different spin parameters to visualize like `ρ`, `T1`, `T2`, among others. 

#### Sequence

There are two options to visualize the sequence already preloaded in RAM: in the time domain or in the k-space. The preloaded sequence is a single-shot EPI.

For visualization of the sequence in the time domain, click on the `Sequence` dropdown and then press the `Sequence (MPS)` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-time-view.png"/></p>
```

For visualization of the sequence in the k-space, click on the `Sequence` dropdown and then press the `k-space` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-kspace-view.png"/></p>
```

You can also display the `Moments` related to the **Sequence** by pressing the `View Moments` and then pressing the buttons for zero, first and second moments.

It is also possible to load **Pulseq** compatible `.seq` sequence files. The **KomaMRI** has some examples stored at `~/.julia/packages/KomaMRI/<id-string>/examples/1.sequences/`. For instance, let's load the `spiral.seq` file and view it the time domain and k-space:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-time-load.png"/></p>
```
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-kspace-load.png"/></p>
```


### Running Simulation
(You can also go to [analog steps using Scripts](programming-workflow.md#Running-Simulation))

Once the inputs are loaded in RAM, it is possible to perform the simulation to get the **Raw Signal**.

#### Simulation Parameters

To visualize the default simulation parameters, click on the `Simulate!` dropdown and then press the `View Options` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-sim-params-view.png"/></p>
```

#### Visualization of the Raw Signal

Press the `Simulate!` button to perform the simulation (this may take a while). Automatically the generated **Raw Signal** should be displayed or you can click on the `Raw Data` dropdown and then press the `View Raw Data` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-rawsignal-view.png"/></p>
```

### Reconstructing Image using MRIReco
(You can also go to [analog steps using Scripts](programming-workflow.md#Reconstructing-Image-using-MRIReco))

Once the **Raw Signal** is loaded in RAM, it is possible to reconstruct the image.

#### Reconstruction Parameters

To visualize the default reconstruction parameters, click on the `Reconstruct!` dropdown and then press the `View Options` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-recon-params-view.png"/></p>
```

#### Visualization of the Image

Press the `Reconstruct!` button to perform the reconstruction (this may take a while).  Automatically the generated **Image** should be displayed or you can click on the he `Reconstruct!` dropdown and then press the `|Image|` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-image-view.png"/></p>
```

### Exporting Results to .mat File
(You can also go to [analog steps using Scripts](programming-workflow.md#Exporting-Results-to-.mat-File))

The user interface has the option to save the results in `.mat` format. Simply click on the `Export to .mat` and you have the alternatives to get data independently or you can press the `All` button to have all the results given by the simulator:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-export-to-mat.png"/></p>
```

So far, and due to limitations of the user interface dependencies, the `.mat` files are saved in the temporal directory of your computer OS, which can be found by typing the `tempdir()` command in the **Julia REPL**.


## Julia Scripts

You should already be familiar with the  [Graphical User Interface](ui-details.md) of **KomaMRI**. However, you can also use this package directly from the **Julia REPL** or write your own Julia scripts. This allows you to unlock the full potential of KomaMRI, enabling you to utilize more of its functionalities and even test your own MRI ideas.

This section demonstrates a basic workflow with **KomaMRI** through writing your own scripts or entering commands directly into the **Julia REPL**. Let's begin.

### Basic Workflow
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

### Loading Simulation Inputs
(You can also go to [analog steps using UI](ui-details.md#Loading-Simulation-Inputs))

The inputs of the simulation are created in the following part of the script: 
```julia
# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
seq = PulseDesigner.EPI_example()
```

#### Scanner

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

#### Phantom

The Phantom struct created in this example represents a slice of a brain. To create it, we use the function `brain_phantom2D`, which is part of the subdependency **KomaMRICore**. While **KomaMRI** provides some phantom examples for experimentation, you may also want to create your custom **Phantom** struct tailored to your specific requirements.

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

You can also visualize the **Phantom** struct using the `plot_phantom_map` function, which is part of the **KomaMRIPlots** subdependency. This function plots the magnitude of a property for each magnetization at a specific spatial position. You can observe properties such as proton density and relaxation times, so feel free to replace the `:ρ` symbol with another property of the phantom in the example below:
```julia-repl
julia> plot_phantom_map(obj, :ρ)
```
```@raw html
<object type="text/html" data="../assets/phantom-rho.html" style="width:100%; height:620px;"></object>
```

To utilize test phantoms included with **KomaMRI**, navigate to the "examples" folder and use the `read_phantom_jemris` function to read a phantom in `.h5` format. The following steps outline how to do this in **Julia**:
```julia-repl
julia> path_koma = dirname(dirname(pathof(KomaMRI)))
julia> path_sphere = joinpath(path_koma, "examples", "2.phantoms", "sphere_chemical_shift.h5")
julia> sphere = read_phantom_jemris(path_sphere)
julia> plot_phantom_map(sphere, :T2)
```
```@raw html
<object type="text/html" data="../assets/phantom-T2-circle.html" style="width:100%; height:620px;"></object>
```

#### Sequence

The **Sequence** struct in the example represents one of the most basic MRI sequences. It excites the object with a 90° RF pulse and then uses EPI gradients to fill the k-space in a "square" manner. While you may want to create your sequences for experiments, you can always use some of the examples already available in **KomaMRI**.

In MRI, the sequence must be carefully designed with precise timing to obtain an image. It includes subcomponents such as gradients, radio-frequency excitation signals, and sample acquisition. For more information on constructing a **Sequence** struct, refer to the [Sequence](sequence.md) section.

You can view general information about a **Sequence** struct by displaying it in the **Julia REPL**:
```julia-repl
julia> seq
Sequence[ τ = 62.846 ms | blocks: 204 | ADC: 101 | GR: 205 | RF: 1 | DEF: 5 ]
```

For more precise timing checks, you can use the `plot_seq` function:
```julia-repl
julia> plot_seq(seq; range=[0 30])
```
```@raw html
<object type="text/html" data="../assets/plot-seq-epi.html" style="width:100%; height:420px;"></object>
```

It is important to consider how the sequence traverses through k-space. The `plot_kspace` function does precisely that:
```julia-repl
julia> plot_kspace(seq)
```
```@raw html
<object type="text/html" data="../assets/kspace-epi.html" style="width:100%; height:420px;"></object>
```

Additionally, there are helpful sequence construction functions within a submodule of **KomaMRI** called **PulseDesigner**. These functions include `RF_hard`, `RF_sinc`, `EPI`, `radial_base` and `spiral_base`. For more details on how to use them, refer to the [API documentation](api.md).

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

### Running Simulation
(You can also go to [analog steps using UI](ui-details.md#Running-Simulation))

The following lines in the example script configure and perform the simulation:
```julia
# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; sim_params)
```

#### Simulation Parameters

To perform simulations, **KomaMRI** requires certain parameters. You can use the default parameters for testing, but you also have the option to customize specific simulation aspects. In the example, we use the `default_sim_params` function to create a dictionary with default simulation parameters. You can view the keys that can be modified by displaying the `sim_params` variable:
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

All of these parameters deserve special attention. We will explain some of the most important ones here. For instance, `"Δt"` and `"Δt_rf"` represent the raster times for the gradients and RFs. `"return_type"` specifies the type of variable returned by the simulator (by default, it returns an object ready for use with **MRIReco** for reconstruction, but you can use the value `"mat"` to return a simple vector). `"gpu"` indicates whether you want to use your GPU device for simulations, and `"precision"` sets the floating-point precision. For more details on how to set these parameters, please refer to the [Simulation Parameters Section](programming-workflow.md#Simulation-Parameters).


#### Raw Signal

The simulation is performed using the `simulate` function, which requires three arguments: a **Scanner** struct, a **Phantom** struct, and a **Sequence** struct. Optionally, you can include the keyword argument `sim_params` if you wish to use custom simulation parameters.

In the example, we can see that the output of the simulation is a special struct:
```julia-repl
julia> typeof(raw)
RawAcquisitionData

julia> raw
RawAcquisitionData[SeqName: epi | 101 Profile(s) of 101×1]
```

You can plot the simulation result with the `plot_signal` function like so:
```julia-repl
julia> plot_signal(raw)
```
```@raw html
<object type="text/html" data="../assets/raw-epi-brain-default.html" style="width:100%; height:420px;"></object>
```


### Reconstructing Image using MRIReco
(You can also go to [analog steps using UI](ui-details.md#Reconstructing-Image-using-MRIReco))

**KomaMRI** does not handle reconstruction; instead, you should utilize the **MRIReco** package to generate an image. For convenience, when you install **KomaMRI**, you also install **MRIReco**, allowing you to access functions from that package. You should pay special attention to the `RawAcquisitionData` and `AcquisitionData` structs, as well as the `reconstruction` function.

In the example below, we define an auxiliary function, `reconstruct_2d_image`, which takes a raw signal struct, `RawAcquisitionData`, as input and returns a 2D Array representing an image. Within this function, we create an `AcquisitionData` struct and set some reconstruction parameters, which serve as inputs for the `reconstruction` function. The latter function is responsible for the image generation process.
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

If you need more information about how to use the `AcquisitionData` and the how to fill the reconstruction parameters, you need to visit the [MRIReco webpage](https://github.com/MagneticResonanceImaging/MRIReco.jl)).

To display the image, you can use the `plot_image` function which is part of the **KomaMRIPlots** subpackage:
```julia-repl
julia> plot_image(image)
```
```@raw html
<center><object type="text/html" data="../assets/image-default-brain.html" style="width:100%; height:620px;"></object></center>
```


### Exporting Results to .mat File
(You can also go to [analog steps using UI](ui-details.md#Exporting-Results-to-.mat-File))

Many people in the MRI community uses MATLAB, probably you are one of them and you want to process the raw signal in the MATLAB environment after simulation is done with **KomaMRI**. Here we show you an example of how to save a `.mat` file with the information of the raw signal thank to the help of the **MAT** package:

Many people in the MRI community use **MATLAB**; you might be one of them and may want to process the **Raw Signal** in the **MATLAB** environment after simulation is completed with **KomaMRI**. Here, we provide an example of how to save a `.mat` file containing the  **Raw Signal** information using the **MAT** package.
```julia
# Use the MAT package
using MAT

# Perform simulation to return an Array type
sim_params["return_type"] = "mat"
raw = simulate(obj, seq, sys; sim_params)

# Save the .mat file in the temp directory
matwrite(joinpath(tempdir(), "koma-raw.mat"), Dict("raw" => raw))
```

Note that we need to simulate to return an array type (not the default `RawAcquisitionData`), and then we utilize the `matwrite` function to save a file named "koma-raw.mat" in your computer's temporary directory. Now, you can navigate to your temporary directory (which you can find by displaying the result of `tempdir()` in the **Julia REPL**) and locate the "koma-raw.mat" file.



## Notebooks

You can use **KomaMRI** with popular programming environments such as **Pluto** and **Jupyter** notebooks. The following sections show how to set up these notebooks and test **KomaMRI** with them.


### Using KomaMRI with Pluto

First, install the **Pluto** module in your **Julia** environment. Remember to press the `]` button to open the **Package Manager Session**:"
```julia-repl
julia>

@(1.9) pkg> add Pluto
```

Afterward, return to the **Julia Session** by pressing the `backspace` button, and then execute the `Pluto.run()` function:
```
julia> using Pluto

julia> Pluto.run()
```

This should automatically open the **Pluto** dashboard in your default web browser:
```@raw html
<center><img src="../assets/pluto-dashboard.png" alt="" style="width:75%;"></center>
```

Next, create a new notebook by clicking on `+ Create a new notebook`:
```@raw html
<center><img src="../assets/pluto-empty-notebook.png" alt="" style="width:75%;"></center>
```

Write and run the following code, which is identical to the [Free Induction Decay example](#Free-Induction-Decay). **Pluto** automatically installs the required modules if they are not present on your system. Additionally, note that we do not use **KomaMRI** directly since we won't be utilizing the `KomaUI()` function. Instead, we rely on the **KomaMRICore** and **KomaMRIPlots** dependencies. To display plots in **Pluto**, ensure that you import the **PlutoPlots** package:"
```@raw html
<embed type="text/html" src="../assets/pluto-fid-example.html" style="height:1024px;width:100%;">
```

### Using KomaMRI with Jupyter

Ensure you have **Jupyter** installed on your computer. Follow [this tutorial](https://test-jupyter.readthedocs.io/en/latest/install.html) for installation using **Anaconda**.

Next, install the **IJulia** module in your **Julia** environment. Remember to press the `]` key to open the **Package Manager Session**:
```julia-repl
julia>

(@v1.9) pkg> add IJulia
```

For this example, make sure to install **KomaMRICore** and **KomaMRIPlots** (we do not use **KomaMRI** directly since we won't be utilizing the `KomaUI()` function):
```julia-repl
(@v1.9) pkg> add KomaMRICore

(@v1.9) pkg> add KomaMRIPlots
```

Next, open **Jupyter**, navigate to a working folder, and create a new notebook by clicking on `New`, then `Julia 1.9.3`."
```@raw html
<center><img src="../assets/jupyter-create-new-notebook.png" alt="" style="width:75%;"></center>
```

A new, empty notebook will appear:
```@raw html
<center><img src="../assets/jupyter-empty-notebook.png" alt="" style="width:75%;"></center>
```

Proceed to write and execute the provided example:

```@raw html
<details><summary>View code</summary>
```
```julia
# Import modules
using KomaMRICore, KomaMRIPlots

# Define sequence
ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF, durRF)

nADC = 8192         # number of acquisition samples
durADC = 250e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay)

seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block

# Plot the sequence
plot_seq(seq; slider=false, height=300)
```
```@raw html
</details>
```

```@raw html
<center><img src="../assets/jupyter-test-notebook.png" alt="" style="width:100%;"></center>
```


This should be sufficient, and now you can start working with **KomaMRI** using **Jupyter** notebooks.

If you encounter the issue of **WebIO** not being detected:
```@raw html
<center><img src="../assets/jupyter-webio-problem.png" alt="" style="width:75%;"></center>
```

Refer to this [troubleshooting guide](https://juliagizmos.github.io/WebIO.jl/stable/troubleshooting/not-detected/) for details. Essentially, you need to install a **WebIO** extension based on your **Jupyter** installation.

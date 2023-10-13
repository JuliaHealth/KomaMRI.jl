# Julia Programming

You should already be familiar with the [Graphical User Interface](ui-details.md) of **KomaMRI.jl**. However, you can always use this package directly from the **Julia REPL** or write your own julia scripts. In this way you can unleash all the power of KomaMRI by using more of its functionalities and you can even test your own MRI ideas.

This section is meant to show you a basic workflow when working with KomaMRI when writing your own scripts or typing down directly on the **Julia REPL**. Let's get started. 

## Basic Workflow

As a general overview, you need to keep in mind the following steps when using KomaMRI:

* Create Inputs: Scanner, Phantom, Sequence
* Perform Simulation
* Perform Reconstruction (actually, this is part of [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl))

Let's replicate these previous steps in a julia script. You will end up with the following, feel free to copy and paste it in the **Julia REPL**:
```julia
# Import the package
using KomaMRI

# Auxiliary function for defining an EPI sequence
function create_epi_seq(sys::Scanner)
    B1 = sys.B1;
    durRF = π/2/(2π*γ*B1)
    EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
    N = 101
    FOV = 23e-2
    EPI = PulseDesigner.EPI(FOV, N, sys)
    TE = 30e-3
    d1 = TE-dur(EPI)/2-dur(EX)
    if d1 > 0 DELAY = Delay(d1) end
    seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI
    seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
    return seq
end

# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
seq = create_epi_seq(sys)

# Define simulation parameters and perform simulation
simParams = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; simParams)

# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqData.traj[1].nodes[:])) #Normalize k-space to -.5 to .5 for NUFFT
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

## Inputs

The inputs of the simulation are created in the following part of the script: 
```julia
# Define scanner, object and sequence
sys = Scanner()
obj = brain_phantom2D()
seq = create_epi_seq(sys)
```

### Scanner

The **Scanner** struct previously created contains default parameters. In your first simulations, you will probably use this default struct without performing any changes to it. You can see all the parameters by simply displaying the struct variable in the **Julia REPL**. The parameters of the Scanner contains hardware limitations like the main magnetic field, the maximum value of the gradients, minimum raster times, among others. These are values that you will be willing to edit in your future custom simulations.  
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

The **Phantom** struct created in this example is an slice of a brain. In order to create it, we use the function **brain\_phantom2D()**, which is part of the subdependency **KomaMRICore**. **KomaMRI** has some phantom examples available for playing around, however maybe you want to create your own **Phantom** struct that fits your needs.

The **Phantom** struct simply contains all the MRI parameters associated with magnetization properties of an object. Those parameters are the positions of the magnetizations, the proton density, the relaxation times, off-resonance, to name some of the most important. You can see all the keys and values of the object in the **julia REPL** like so:
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
As you can see, the attributes of the **Phantom** struct are vectors of the object properties in which each element contains a value associated to a single magnetization.

You can also display the **Phantom** struct with the help of the function **plot\_phantom\_map()**, which is part of the **KomaMRIPlots** subdependency. The function to plot the phantom displays the magnitude of a property for each magnetization at a certain space position. You can see properties like the proton density or the relaxation times, so feel free to change the **:ρ** symbol by other property of the phantom in the example bellow.
```julia-repl
julia> plot_phantom_map(obj, :ρ)
```
```@raw html
<p align="center"><img width="100%" src="../assets/phantom-rho.svg"/></p>
```

For using some test phantoms installed with **KomaMRI**, you need to navigate in the "examples" folder and then use the function **read\_phantom\_jemris()** to read a phantom in **.h5** format. Here are shown the steps of how is done in **Julia**:
```julia-repl
julia> path_koma = dirname(dirname(pathof(KomaMRI)))
julia> path_sphere = joinpath(path_koma, "examples", "2.phantoms", "sphere_chemical_shift.h5")
julia> sphere = read_phantom_jemris(path_sphere)
julia> plot_phantom_map(sphere)
```
```@raw html
<p align="center"><img width="100%" src="../assets/phantom-T2-circle.svg"/></p>
```

### Sequence

The **Sequence** struct in the example represents one of the most basic MRI sequences. It excites the object with a 90° RF pulse and then uses EPI gradients to fill the k-space in an ''square'' manner. For sure you want to create your own sequences in your experiments, however you always can use some examples already available in **KomaMRI**.

In MRI, the sequence must be designed carefully and with precise times in order to obtain an image. It contains other sub components like gradients, radio-frequency excitation signals and acquisition of samples. More about how to construct a **Sequence** struct is explained in the section [Sequence](sequence.md).

You can display some general information about a **Sequence** struct by simply displaying it on the **Julia REPL**:
```julia-repl
julia> seq
Sequence[ τ = 62.846 ms | blocks: 204 | ADC: 101 | GR: 205 | RF: 1 | DEF: 5 ]
```

Or even more helpful for checking precise timings, you can use the **plot\_seq()** function.
```julia-repl
julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="100%" src="../assets/plot-seq-epi.svg"/></p>
```

A helpful information about the sequence is to keep in mind the way how it travels through the k-space. The function **get\_kspace()** does exactly that:
```julia-repl
julia> plot_kspace(seq)
```
```@raw html
<p align="center"><img width="100%" src="../assets/kspace-epi.svg"/></p>
```

Additionally, there are some helpful functions for constructing sequences which are part of a submodule of **KomaMRI** called **PulseDesigner**. Functions like **RF\_hard()**, **RF\_sinc()**, **EPI()**, **radial\_base()** and **spiral\_base()**. Check the [API documentation](api.md) for more details about how to use them.

**KomaMRI** also offers compatibility with [Pulseq](https://pulseq.github.io/). The package installation also has some **.seq** files (**Pulseq** format) which can be read and use them as a **Sequence** struct. Here is shown how you can read a spiral **Pulseq** file which is stored in the "examples" folder of **KomaMRI**:
```julia-repl
julia> path_koma = dirname(dirname(pathof(KomaMRI)))
julia> path_spiral = joinpath(path_koma, "examples", "1.sequences", "spiral.seq")
julia> spiral = read_seq(path_spiral)
julia> plot_seq(spiral)
julia> plot_kspace(spiral)
```
```@raw html
<p align="center"><img width="100%" src="../assets/seq-spiral-pulseq.svg"/></p>
<p align="center"><img width="100%" src="../assets/seq-spiral-pulseq-kspace.svg"/></p>
```

## Simulation

The lines in the example script where the simulation is configured and performed are the following: 
```julia
# Define simulation parameters and perform simulation
simParams = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; simParams)
```

### Simulation Parameters

In order to simulate, **KomaMRI** requires some parameters. You can always use the default parameters for testing purposes, however you have the option of customize certain simulation aspects. In the example, we use the function **default\_sim\_params()** to create a dictionary with default simulation parameters. You can see the keys that can be modified by displaying the "simParams" variable:
```julia-repl
julia> simParams
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

All of these parameters deserves special attention, however we are going to explain some of the most important here. For instance, "Δt" and "Δt_rf" are the raster times for the gradients and RFs, "return_type" is the type of variable that is returned by the simulator (by default it returns an object that is ready to be used with **MRIReco** for reconstruction, you can use the value "mat" to return a simple vector type), "gpu" indicates if you want yo use your GPU device to perform simulations, "precision" is for floating point precision.

### Raw Signal

The simulation is performed with the help of the **simulate()** function. As you can notice, it requires 3 arguments: a **Scanner** struct, a **Phantom** struct and a **Sequence** struct. Optionally, you can add the keyword argument "simParams" if you want to perform a simulation with default simulation parameters. 

In the example, the output of the simulation is a special struct:
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
<p align="center"><img width="100%" src="../assets/raw-epi-brain-default.svg"/></p>
```

Many people in the MRI community uses MATLAB, probably you are one of them and you want to process the raw signal in the MATLAB environment. Here we show you an example of how to save a ".mat" file with the information of raw signal thank to the help of the **MAT** package:
```julia
# Use the MAT package
using MAT

# Perform simulation to return an Array type
simParams["return_type"] = "mat"
raw = simulate(obj, seq, sys; simParams)

# Save the .mat file in the temp directory
matwrite(joinpath(tempdir(), "koma-raw.mat"), Dict("raw" => raw))
```

Note that we need to simulate again to return an array type and then we use the **matwrite()** function to save a file called "koma-raw.mat" in the temporal directory of the your computer OS. Now you can go to your temporal directory (you can check that by displaying the **tempdir()** result in the **Julia REPL**) and see the "koma-raw.mat" file.


## Reconstruction

**KomaMRI** doesn't perform reconstruction, instead you need to use the **MRIReco** package to obtain an image. For simplicity, when you install **KomaMRI**, you also install **MRIReco**, so you can use functions from that package too. In particular, you will need to pay attention to the structs **RawAcquisitionData**, **AcquisitionData** and to the function **reconstruction()**.

In the example bellow, we create a auxiliary function **reconstruct\_2d\_image()** which receives as an input a raw signal struct **RawAcquisitionData** and returns a 2D Array which represents a 2D image. Inside this function, we create a **AcquisitionData** struct and also we fill some reconstruction parameters which are inputs of the **reconstruction()** function. This last function is in charge of performing the process of getting the image. 
```julia
# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqData.traj[1].nodes[:])) #Normalize k-space to -.5 to .5 for NUFFT
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
<p align="center"><img width="100%" src="../assets/image-default-brain.svg"/></p>
```

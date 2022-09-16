# Simulation Examples

In order to follow these examples, make sure you have installed Julia in your machine and the `KomaMRI` package. These are step-by-step examples, so you can simply copy-and-paste the next lines of code in the Julia REPL to follow along.

Before starting any example, don't forget to Include the `KomaMRI` package:
```
using KomaMRI
```

## Brain Example

This example is the same as the one used by default in the [KomaMRI User Interface](getting-started.md#UI-Example) but instead of using the **KomaUI()** function we issue commands in the Julia REPL with functions from the [KomaMRI API](api.md). The procedures in this example can be extended to other custom application cases by simply changing the inputs (**Scanner**, **Sequence** and **Phantom**).
```@setup 0
using KomaMRI
```

### Preparation: Inputs

Define the default **Scanner** struct. Note that it simply contains some parameters that defines a scanner device:
```@example 0
sys = Scanner()
```

Define the **Sequence** struct. In this case we use a one-shot 90-degree excitation and a EPI sequence. It is also a good idea to inspect how the signal looks like in the time domain and how it covers the kspace:
```@setup 0
# Define some parameters
FOV, N = 23e-2, 101         # field-of-view and number-of-pixels in x,y
durRF = π/2/(2π*γ*sys.B1)   # duration of a flat RF for 90-degree hard excitation pulse

# Define the sequence
ex = PulseDesigner.RF_hard(sys.B1, durRF, sys)  # the 90-degree excitation
epi = PulseDesigner.EPI(FOV, N, sys)            # the epi sequence
seq = ex + epi                                  # the final sequence

# Plot the sequence in time and in the kspace
savefig(plot_seq(seq; slider=true, height=300), "assets/0-seq.html")
savefig(plot_kspace(seq; height=700), "assets/0-kspace.html")
```
```julia
# Define some parameters
FOV, N = 23e-2, 101         # field-of-view and number-of-pixels in x,y
durRF = π/2/(2π*γ*sys.B1)   # duration of a flat RF for 90-degree hard excitation pulse

# Define the sequence
ex = PulseDesigner.RF_hard(sys.B1, durRF, sys)  # the 90-degree excitation
epi = PulseDesigner.EPI(FOV, N, sys)            # the epi sequence
seq = ex + epi                                  # the final sequence

# Plot the sequence in time and in the kspace
plot_seq(seq)
plot_kspace(seq)
```
```@raw html
<object type="text/html" data="../assets/0-seq.html" style="width:100%; height:330px;"></object>
```
```@raw html
<object type="text/html" data="../assets/0-kspace.html" style="width:100%; height:730px;"></object>
```

Define the **Phantom** struct. This is an example of a 2D brain. We can even visualize the parameters of the spins in an image, for example here we plot the densities of every spin in the phantom:
```@setup 0
# Define the phantom
obj = brain_phantom2D()     # an example of a 2D brain

# Visualize the densities of the spins
savefig(plot_phantom_map(obj, :ρ; height=700), "assets/0-obj.html")
```
```julia
# Define the phantom
obj = brain_phantom2D()     # an example of a 2D brain

# Visualize the densities of the spins
plot_phantom_map(obj, :ρ)
```
```@raw html
<object type="text/html" data="../assets/0-obj.html" style="width:100%; height:730px;"></object>
```

### Simulation: Raw Signal Output

Simulate and get the raw signal in ISMRMRD format:
```@setup 0
# Simulate and get the raw signal in ISMRMRD format
signal = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys)

# Plot the raw signal
savefig(plot_signal(ismrmrd; height=300), "assets/0-ismrmrd.html")
```
```julia
# Simulate and get the raw signal in ISMRMRD format
signal = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys)

# Plot the raw signal
plot_signal(ismrmrd)
```
```@raw html
<object type="text/html" data="../assets/0-ismrmrd.html" style="width:100%; height:330px;"></object>
```

### Reconstruction: Image Output

To reconstruct the image we use the [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl) package. Here it is an example of how to use the raw signal in ISMRMRD format with **MRIReco.jl** to finally visualize the image:
```@setup 0
# Get the acquisition data
Nx, Ny = ismrmrd.params["reconSize"][1:2]
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true)
acq = AcquisitionData(ismrmrd)

# Reconstruct and get an slice of the magnitude of the image
recon = reconstruction(acq, params)
image  = reshape(recon.data, Nx, Ny, :)
slice_abs = abs.(image[:, :, 1])

# Plot an slice of the image
savefig(plot_image(slice_abs; height=700), "assets/0-slice_abs.html")
```
```julia
# Get the acquisition data
Nx, Ny = ismrmrd.params["reconSize"][1:2]
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true)
acq = AcquisitionData(ismrmrd)

# Reconstruct and get an slice of the magnitude of the image
recon = reconstruction(acq, params)
image  = reshape(recon.data, Nx, Ny, :)
slice_abs = abs.(image[:, :, 1])

# Plot an slice of the image
plot_image(slice_abs)
```
```@raw html
<object type="text/html" data="../assets/0-slice_abs.html" style="width:100%; height:730px;"></object>
```


## 1-Spin Example

This example simulates a phantom that contain only 1 spin, so it is a very theoretical example. The idea is to visualize the dynamics of an spin when an RF pulse is applied. Two cases are considered: the simple exponential free decay and when the spin is affected by the off-resonance phenomena: 

```@setup 1
using KomaMRI
```

### Sequence Definition

Let's define the default scanner struct:
```@setup 1
sys = Scanner()
```
```julia
sys = Scanner()
```

In this example we define a sequence with 3 blocks. The first block contains a 90° RF pulse, the second block is a small delay and the third one is the adquisition (for more information about how a sequence is constructed refer to the [Sequence Structure](useful-information.md#Sequence-Structure)):
```@setup 1
# Define the 90-degree hard excitation pulse parameters
ampRF = sys.B1                      # amplitude of the RF pulse
durRF = π / 2 / (2π * γ * ampRF)    # duration of the RF pulse (γ is globally defined by KomaMRI)

# Define a delay parameter between the excitation and the acquisition
delayTime = .001    # a small delay

# Define parameters for the acquisition
nADC = 8192         # number of acquisition samples
durADC = .5         # duration of the acquisition

# Define sequence the struct
exc = Sequence([Grad(0,0); Grad(0,0); Grad(0,0);;], [RF(ampRF,durRF);;])
dly = Delay(delayTime)
adq = Sequence([Grad(0,0); Grad(0,0); Grad(0,0);;], [RF(0,0);;], [ADC(nADC, durADC)])
seq = exc + dly + adq

# Plot the sequence over time
savefig(plot_seq(seq; slider=true, height=300), "assets/1-seq.html")
```
```julia
# Define the 90-degree hard excitation pulse parameters
ampRF = sys.B1                      # amplitude of the RF pulse
durRF = π / 2 / (2π * γ * ampRF)    # duration of the RF pulse (γ is globally defined by KomaMRI)

# Define a delay parameter between the excitation and the acquisition
delayTime = durRF   # a small delay

# Define parameters for the acquisition
nADC = 8192         # number of acquisition samples
durADC = .5         # duration of the acquisition

# Define sequence the struct
exc = Sequence([Grad(0,0); Grad(0,0); Grad(0,0);;], [RF(ampRF,durRF);;])
dly = Delay(delayTime)
adq = Sequence([Grad(0,0); Grad(0,0); Grad(0,0);;], [RF(0,0);;], [ADC(nADC, durADC)])
seq = exc + dly + adq

# Plot the sequence over time
plot_seq(seq; slider=true)
```
```@raw html
<object type="text/html" data="../assets/1-seq.html" style="width:100%; height:330px;"></object>
```

### Free-Decay Case

Define the 1-spin phantom struct with T1 and T2 parameters:
```@setup 1
# Define the 1-spin phantom struct with T1 and T2 parameters
T1, T2 = 1, durADC / 5
obj = Phantom(x=[0], T1=[T1], T2=[T2])
```
```julia
# Define the 1-spin phantom struct with T1 and T2 parameters
T1, T2 = 1, durADC / 5
obj = Phantom(x=[0], T1=[T1], T2=[T2])
```

Perform the simulation and visualize the free-decay behavior:
```@setup 1
# Simulate and get the output raw signal in ismrmrd format
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)

# Check that the raw signal is the expected free-decay behavior
savefig(plot_signal(ismrmrd; height=300), "assets/10-ismrmrd.html")
```
```julia
# Simulate and get the output raw signal in ismrmrd format
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)

# Check that the raw signal is the expected free-decay behavior
plot_signal(ismrmrd)
```
```@raw html
<object type="text/html" data="../assets/10-ismrmrd.html" style="width:100%; height:330px;"></object>
```

### Off-Resonance Case

Define the 1-spin phantom struct with the additional off-resonance parameter:
```@setup 1
Δw = 220 * 2π
obj = Phantom(x=[0], T1=[T1], T2=[T2], Δw=[Δw])
```
```julia
Δw = 220 * 2π
obj = Phantom(x=[0], T1=[T1], T2=[T2], Δw=[Δw])
```

Perform the simulation and visualize the off-resonance oscillation behavior:
```@setup 1
# Simulate and get the output raw signal in ismrmrd format
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)

# Check that the raw signal is the expected off-resonance oscillation behavior
savefig(plot_signal(ismrmrd; height=300), "assets/11-ismrmrd.html")
```
```julia
# Simulate and get the output raw signal in ismrmrd format
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)

# Check that the raw signal is the expected off-resonance oscillation behavior
plot_signal(ismrmrd)
```
```@raw html
<object type="text/html" data="../assets/11-ismrmrd.html" style="width:100%; height:330px;"></object>
```


## Shifted Sphere

This example shows the simulation and reconstruction for a sphere phantom with a chemical shift. The idea is to include standard files to define the **Sequence** and **Phantom** inputs . The files can be found in the same paths of the [github repo](https://github.com/cncastillo/KomaMRI.jl).
```@setup 2
using KomaMRI
```

### Preparation: Inputs

Define the inputs of the simulation (the scanner, the sequence and the phantom):
```@setup 2
# Define the scanner, sequence and phantom
sys = Scanner()
seq = read_seq("../../examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq")
obj = read_phantom_jemris("../../examples/2.phantoms/sphere_chemical_shift.h5")

# Visualize the inputs
display(sys)
savefig(plot_seq(seq; slider=true, height=300), "assets/2-seq.html")
savefig(plot_kspace(seq; height=700), "assets/2-kspace.html")
savefig(plot_phantom_map(obj, :T2; height=700), "assets/2-obj.html")
```
```julia
# Define the scanner, sequence and phantom
sys = Scanner()
seq = read_seq("examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq")
obj = read_phantom_jemris("examples/2.phantoms/sphere_chemical_shift.h5")

# Visualize the inputs
display(sys)
plot_seq(seq; slider=true)
plot_kspace(seq)
plot_phantom_map(obj, :T2)
```
```@raw html
<object type="text/html" data="../assets/2-seq.html" style="width:100%; height:330px;"></object>
```
```@raw html
<object type="text/html" data="../assets/2-kspace.html" style="width:100%; height:730px;"></object>
```
```@raw html
<object type="text/html" data="../assets/2-obj.html" style="width:100%; height:730px;"></object>
```

### Simulation: Raw Signal Output

Simulate and get the raw signal in ISMRMRD format:
```@setup 2
# Simulate and get the raw signal in ISMRMRD format
signal = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys)

# Plot the raw signal
savefig(plot_signal(ismrmrd; height=300), "assets/2-ismrmrd.html")
```
```julia
# Simulate and get the raw signal in ISMRMRD format
signal = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys)

# Plot the raw signal
plot_signal(ismrmrd)
```
```@raw html
<object type="text/html" data="../assets/2-ismrmrd.html" style="width:100%; height:330px;"></object>
```


### Reconstruction: Image Output

To reconstruct the image we use the [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl) package:
```@setup 2
# Get the acquisition data
Nx, Ny = ismrmrd.params["reconSize"][1:2]
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true)
acq = AcquisitionData(ismrmrd)

# Reconstruct and get an slice of the magnitude of the image
recon = reconstruction(acq, params)
image  = reshape(recon.data, Nx, Ny, :)
slice_abs = abs.(image[:, :, 1])

# Plot an slice of the image
savefig(plot_image(slice_abs; height=700), "assets/2-slice_abs.html")
```
```julia
# Get the acquisition data
Nx, Ny = ismrmrd.params["reconSize"][1:2]
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true)
acq = AcquisitionData(ismrmrd)

# Reconstruct and get an slice of the magnitude of the image
recon = reconstruction(acq, params)
image  = reshape(recon.data, Nx, Ny, :)
slice_abs = abs.(image[:, :, 1])

# Plot an slice of the image
plot_image(slice_abs)
```
```@raw html
<object type="text/html" data="../assets/2-slice_abs.html" style="width:100%; height:730px;"></object>
```

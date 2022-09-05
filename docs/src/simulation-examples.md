# Simulation Examples

In order to follow these examples, make sure you have installed Julia in your machine and the `KomaMRI` package. For more advanced examples, you also need to download or clone this repo in your local machine since we are going to import some example files to generate phantom and sequence objects. These are step-by-step examples, so you can simply copy-and-paste the next lines of code in the Julia REPL to follow along.

## 1-Spin Example

Include the `KomaMRI` package:
```@example 0
using KomaMRI
```

Let's define the default scanner struct. This is an input for the simulator and contains some parameters that defines an scanner device:
```@example 0
sys = Scanner()
```

Other input for the simulator is the sequence struct. A sequence can be though as and ordered concatenation of blocks over time. Every block is composed by the x-y-z gradients, an RF pulse and the acquisition of the samples, which are also defined by the structs GR, RF and ADC respectively. For gradients and RF pulses structs, it's necessary to define amplitude and duration parameters. For the adquisition, it must be defined the number of samples and its duration. Note that sequence blocks can have different time durations.

In this example we define a sequence with 3 blocks. The first block contains a 90° RF pulse, the second block contains the acquisition and the third block is just an additional block with small gradient amplitudes and a small duration. Note that every sequence block must define the GR, RF and ADC structs, but when some of them are not necessary to consider we simply define a zero first argument (zero amplitud for GR and RF or zero number of samples for ADC):
```@setup 0
# Define the 90-degree hard excitation pulse parameters
ampRF = sys.B1                      # amplitude of the RF pulse
durRF = π / 2 / (2π * γ * ampRF)    # duration of the RF pulse (γ is globally defined by KomaMRI)

# Define the 90-degree hard excitation pulse parameters
nADC = 128          # number of acquisition samples
durADC = .5         # duration of the acquisition
delayADC = .001     # a small delay

# Define the parameters for the last block of the sequence
ampGR = .1 * ampRF      # amplitude of the gradient in the last block
durGR = .1 * durRF      # duration of the gradient in the last block

# Define sequence struct
matGRs  = [Grad(0, durRF)    Grad(0, durADC)               Grad(ampGR, durGR);
           Grad(0, durRF)    Grad(0, durADC)               Grad(ampGR, durGR);
           Grad(0, durRF)    Grad(0, durADC)               Grad(ampGR, durGR)]
matRFs  = [RF(ampRF, durRF)  RF(0, durADC)                 RF(0, durGR)]
vecADCs = [ADC(0, durRF);    ADC(nADC, durADC, delayADC);  ADC(0, durGR)]
seq = Sequence(matGRs, matRFs, vecADCs)
```
```
# Define the 90-degree hard excitation pulse parameters
ampRF = sys.B1                      # amplitude of the RF pulse
durRF = π / 2 / (2π * γ * ampRF)    # duration of the RF pulse (γ is globally defined by KomaMRI)

# Define the 90-degree hard excitation pulse parameters
nADC = 1000          # number of acquisition samples
durADC = .5         # duration of the acquisition
delayADC = .001     # a small delay

# Define the parameters for the last block of the sequence
ampGR = .1 * ampRF      # amplitude of the gradient in the last block
durGR = .1 * durRF      # duration of the gradient in the last block

# Define sequence struct
matGRs  = [Grad(0, durRF)    Grad(0, durADC)               Grad(ampGR, durGR);
           Grad(0, durRF)    Grad(0, durADC)               Grad(ampGR, durGR);
           Grad(0, durRF)    Grad(0, durADC)               Grad(ampGR, durGR)]
matRFs  = [RF(ampRF, durRF)  RF(0, durADC)                 RF(0, durGR)]
vecADCs = [ADC(0, durRF);    ADC(nADC, durADC, delayADC);  ADC(0, durGR)]
seq = Sequence(matGRs, matRFs, vecADCs)
```

We can visualize the created input sequence:
```@setup 0
p = plot_seq(seq; slider=true, height=300)
savefig(p, "assets/0-seq.html")
```
```
plot_seq(seq; slider=true, height=300)
```
```@raw html
<object type="text/html" data="../assets/0-seq.html" style="width:100%; height:330px;"></object>
```

The last remaining input to define is the phantom. In this simple example we are going to define a 1-spin phantom with T1 and T2 decay parameters, and additionally we are going to see the effect when we consider the off-resonance effect. Before going further, we need to define the voxel grid Nx-Ny-Nz parameters, since this is just 1-spin example, we simply assign 1 voxel.
```@setup 0
seq.DEF["Nx"], seq.DEF["Ny"], seq.DEF["Nz"] = 1, 1, 1
```
```
seq.DEF["Nx"], seq.DEF["Ny"], seq.DEF["Nz"] = 1, 1, 1
```

### Free Decay

Define the 1-spin phantom struct with T1 and T2 parameters:
```@setup 0
T1, T2 = 1, durADC / 5
obj = Phantom(x=[0], T1=[T1], T2=[T2])
```
```
T1, T2 = 1, durADC / 5
obj = Phantom(x=[0], T1=[T1], T2=[T2])
```

Perform the simulation and transform the output raw signal to ismrmrd format:
```@setup 0
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)
```
```
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)
```

The output raw signal is the expected free decay behavior:
```@setup 0
p = plot_signal(ismrmrd; height=300)
savefig(p, "assets/00-ismrmrd.html")
```
```
plot_signal(ismrmrd; height=300)
```
```@raw html
<object type="text/html" data="../assets/00-ismrmrd.html" style="width:100%; height:330px;"></object>
```

### Off-Resonance

Define the 1-spin phantom struct with the additional off-resonance parameter:
```@setup 0
Δw = 220 * 2π
obj = Phantom(x=[0], T1=[T1], T2=[T2], Δw=[Δw])
```
```
Δw = 220 * 2π
obj = Phantom(x=[0], T1=[T1], T2=[T2], Δw=[Δw])
```

Perform the simulation and transform the output raw signal to ismrmrd format:
```@setup 0
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)
```
```
sig = simulate(obj, seq, sys)
ismrmrd = rawSignalToISMRMRD([sig;;], seq; phantom=obj, sys=sys)
```

The output raw signal is the expected off-resonance oscillation behavior:
```@setup 0
p = plot_signal(ismrmrd; height=300)
savefig(p, "assets/01-ismrmrd.html")
```
```
plot_signal(ismrmrd; height=300)
```
```@raw html
<object type="text/html" data="../assets/01-ismrmrd.html" style="width:100%; height:330px;"></object>
```

## `sphere_chemical_shift`

### Simulation

Include the `KomaMRI` package:
```@example 1
using KomaMRI
```

Define the inputs of the simulation (the phantom, the sequence and the scanner):
```@setup 1
obj = read_phantom_jemris("../../examples/2.phantoms/sphere_chemical_shift.h5");
seq = read_seq("../../examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq");
sys = Scanner();
```
```
obj = read_phantom_jemris("examples/2.phantoms/sphere_chemical_shift.h5");
seq = read_seq("examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq");
sys = Scanner();
```

Simulate and get the raw signal in ISMRMRD format:
```@setup 1
signal = simulate(obj, seq, sys);
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys);
```
```
signal = simulate(obj, seq, sys);
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys);
```

Plot the input sequence:
```@setup 1
p = plot_seq(seq; slider=true, height=300);
savefig(p, "assets/1-seq.html");
```
```
plot_seq(seq; slider=true, height=300)
```
```@raw html
<object type="text/html" data="../assets/1-seq.html" style="width:100%; height:330px;"></object>
```

Plot the k-space:
```@setup 1
p = plot_kspace(seq; height=700);
savefig(p, "assets/1-kspace.html");
```
```
plot_kspace(seq; height=700)
```
```@raw html
<object type="text/html" data="../assets/1-kspace.html" style="width:100%; height:730px;"></object>
```

Plot the output raw signal:
```@setup 1
p = plot_signal(ismrmrd; height=300);
savefig(p, "assets/1-ismrmrd.html");
```
```
plot_signal(ismrmrd; height=300)
```
```@raw html
<object type="text/html" data="../assets/1-ismrmrd.html" style="width:100%; height:330px;"></object>
```

### Reconstruction

We need to define some variables first:
```@setup 1
Nx, Ny = ismrmrd.params["reconSize"][1:2];
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true);
acq = AcquisitionData(ismrmrd);
```
```
Nx, Ny = ismrmrd.params["reconSize"][1:2];
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true);
acq = AcquisitionData(ismrmrd);
```

Then we can reconstruct the image an get the magnitude for an slice of it (in general the image is a 3 dimensional array):
```@setup 1
recon = reconstruction(acq, params);
image  = reshape(recon.data, Nx, Ny, :);
slice_abs = abs.(image[:, :, 1]);
```
```
recon = reconstruction(acq, params);
image  = reshape(recon.data, Nx, Ny, :);
slice_abs = abs.(image[:, :, 1]);
```

Plot the magnitude of the image slice:
```@setup 1
p = plot_image(slice_abs; height=700);
savefig(p, "assets/1-slice_abs.html");
```
```
plot_image(slice_abs; height=700)
```
```@raw html
<object type="text/html" data="../assets/1-slice_abs.html" style="width:100%; height:730px;"></object>
```


## `brain_susceptibility`

### Simulation

Include the `KomaMRI` package:
```@example 2
using KomaMRI
```

Define the inputs of the simulation (the phantom, the sequence and the scanner):
```@setup 2
obj = read_phantom_jemris("../../examples/2.phantoms/brain_susceptibility.h5");
seq = read_seq("../../examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq");
sys = Scanner();
```
```
obj = read_phantom_jemris("examples/2.phantoms/brain_susceptibility.h5");
seq = read_seq("examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq");
sys = Scanner();
```

Simulate and get the raw signal in ISMRMRD format:
```@setup 2
signal = simulate(obj, seq, sys);
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys);
```
```
signal = simulate(obj, seq, sys);
ismrmrd = rawSignalToISMRMRD([signal;;], seq; phantom=obj, sys=sys);
```

Plot the input sequence:
```@setup 2
p = plot_seq(seq; slider=true, height=300);
savefig(p, "assets/2-seq.html");
```
```
plot_seq(seq; slider=true, height=300)
```
```@raw html
<object type="text/html" data="../assets/2-seq.html" style="width:100%; height:330px;"></object>
```

Plot the k-space:
```@setup 2
p = plot_kspace(seq; height=700);
savefig(p, "assets/2-kspace.html");
```
```
plot_kspace(seq; height=700)
```
```@raw html
<object type="text/html" data="../assets/2-kspace.html" style="width:100%; height:730px;"></object>
```

Plot the output raw signal:
```@setup 2
p = plot_signal(ismrmrd; height=300);
savefig(p, "assets/2-ismrmrd.html");
```
```
plot_signal(ismrmrd; height=300)
```
```@raw html
<object type="text/html" data="../assets/2-ismrmrd.html" style="width:100%; height:330px;"></object>
```

### Reconstruction

We need to define some variables first:
```@setup 2
Nx, Ny = ismrmrd.params["reconSize"][1:2];
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true);
acq = AcquisitionData(ismrmrd);
```
```
Nx, Ny = ismrmrd.params["reconSize"][1:2];
params = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny), :densityWeighting=>true);
acq = AcquisitionData(ismrmrd);
```

Then we can reconstruct the image an get the magnitude of an slice of it (in general the image is a 3 dimensional array):
```@setup 2
recon = reconstruction(acq, params);
image  = reshape(recon.data, Nx, Ny, :);
slice_abs = abs.(image[:, :, 1]);
```
```
recon = reconstruction(acq, params);
image  = reshape(recon.data, Nx, Ny, :);
slice_abs = abs.(image[:, :, 1]);
```

Plot the magnitude of the image slice:
```@setup 2
p = plot_image(slice_abs; height=700);
savefig(p, "assets/2-slice_abs.html");
```
```
plot_image(slice_abs; height=700)
```
```@raw html
<object type="text/html" data="../assets/2-slice_abs.html" style="width:100%; height:730px;"></object>
```

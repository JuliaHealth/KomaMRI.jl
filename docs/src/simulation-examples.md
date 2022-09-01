# Simulation Examples

In order to follow these examples, make sure you have installed Julia in your machine and the `KomaMRI` package. Also, you need to download or clone this repo in your local machine since we are going to import some example files to generate phantom and sequence objects.

## 1-spin Example


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
plot_seq(seq; slider=true, height=300);
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
plot_kspace(seq; height=700);
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
plot_signal(ismrmrd; height=300);
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
plot_image(slice_abs; height=700);
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
plot_seq(seq; slider=true, height=300);
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
plot_kspace(seq; height=700);
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
plot_signal(ismrmrd; height=300);
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
plot_image(slice_abs; height=700);
```
```@raw html
<object type="text/html" data="../assets/2-slice_abs.html" style="width:100%; height:730px;"></object>
```

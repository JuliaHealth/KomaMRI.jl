# Graphical User Interface

This section is meant to explain some details about how to use the user interface of the **KomaMRI.jl** package and the processes that happen internally while interacting with it.

## Launch the UI

Open the Julia REPL and issue the following commands to include the **KomaMRI.jl** package and launch the user interface:
```julia-repl
julia> using KomaMRI

julia> KomaUI()
```
```@raw html
<p align="center"><img width="90%" src="../assets/gui-dashboard.png"/></p>
```

## Inputs

The user interface has already preloaded some inputs (stored in RAM). In particular, it has predefined the **Scanner**, the **Phantom** and the **Sequence** structs. In the following subsections, we will show how to visualize these inputs.

### Scanner

So far, it is not possible to see the **Scanner** struct in the user interface. However, the preloaded **Scanner** struct is the default one, so it is possible to know its attributes by creating a new default **Scanner** struct in the Julia REPL like so:
```julia-repl
julia> sys = Scanner()
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

To see the phantom already stored in RAM, simply click on the `Phantom` dropdown an then press the `View Phantom` button. The preloaded phantom is a slice of a brain:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-phantom-view.png"/></p>
```

It is also possible to load `.h5` phantom files. The **KomaMRI.jl** has some examples stored at `~/.julia/packages/KomaMRI/<id-string>/examples/2.phantoms/`. For instance, let's load the `sphere_chemical_shift.h5` file:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-phantom-load.png"/></p>
```

Note that you can select different spin parameters to visualize like ``\rho``, ``T_1``, ``T_2``, among others. 

### Sequence

There are two options to visualize the sequence already preloaded in RAM: in the time domain or in the k-space. The preloaded sequence is a single-shot EPI.

For visualization of the sequence in the time domain, click on the `Sequence` dropdown and then press the `Sequence (MPS)` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-time-view.png"/></p>
```

For visualization of the sequence in the k-space, click on the `Sequence` dropdown and then press the `k-space` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-kspace-view.png"/></p>
```

It is also possible to load `.seq` sequence files. The **KomaMRI.jl** has some examples stored at `~/.julia/packages/KomaMRI/<id-string>/examples/1.sequences/`. For instance, let's load the `spiral.seq` file and view it the time domain and k-space:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-time-load.png"/></p>
```
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-kspace-load.png"/></p>
```


## Simulation

Once the inputs are loaded in RAM, it is possible to perform the simulation to get the **Raw Signal**.

### Simulation Parameters

To visualize the default simulation parameters, click on the `Simulate!` dropdown and then press the `View Options` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-sim-params-view.png"/></p>
```

### Visualization of the Raw Signal

Press the `Simulate!` button to perform the simulation (this may take a while). Then, to view the generated **Raw Signal**, click on the `Raw Data` dropdown and then press the `View Raw Data` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-rawsignal-view.png"/></p>
```

## Reconstruction

Once the **Raw Signal** is loaded in RAM, it is possible to reconstruct the image.

### Reconstruction Parameters

To visualize the default reconstruction parameters, click on the `Reconstruct!` dropdown and then press the `View Options` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-recon-params-view.png"/></p>
```

### Visualization of the Image

Press the `Reconstruct!` button to perform the reconstruction (this may take a while). Then, to view the generated **Image**, click on the he `Reconstruct!` dropdown and then press the `|Image|` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-image-view.png"/></p>
```
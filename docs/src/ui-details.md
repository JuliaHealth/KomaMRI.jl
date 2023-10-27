# Graphical User Interface

This section is meant to explain some details about how to use the user interface of the **KomaMRI** package and the processes that happen internally while interacting with it.

## Basic Workflow
(You can also go to [analog steps using Scripts](programming-workflow.md#Basic-Workflow))

As a general overview, you need to keep in mind the following steps when using KomaMRI:

* Loading Simulation Inputs: Scanner, Phantom, Sequence
* Running Simulation
* Reconstructing Image using MRIReco

In the following subsections, we will cover all the aforementioned steps. But first, open the Julia REPL and issue the following commands to include the **KomaMRI** package and launch the user interface:
```julia-repl
julia> using KomaMRI

julia> KomaUI()
```
```@raw html
<p align="center"><img width="90%" src="../assets/gui-dashboard.png"/></p>
```

## Loading Simulation Inputs
(You can also go to [analog steps using Scripts](programming-workflow.md#Loading-Simulation-Inputs))

The user interface has already preloaded some inputs (stored in RAM). In particular, it has predefined the **Scanner**, the **Phantom** and the **Sequence** structs. In the following subsections, we will show how to visualize these inputs.

### Scanner

You can visualize the preloaded **Scanner** struct by clicking on the `Scanner` dropdown and then pressing the `View Scanner` button. The **Scanner** struct simply contains information related to the hardware being use like the magnitude of the main magnetic field:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-scanner-view.png"/></p>
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

You can also display the `Moments` related to the **Sequence** by pressing the `View Moments` and then pressing the buttons for zero, first and second moments.

It is also possible to load `.seq` sequence files. The **KomaMRI** has some examples stored at `~/.julia/packages/KomaMRI/<id-string>/examples/1.sequences/`. For instance, let's load the `spiral.seq` file and view it the time domain and k-space:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-time-load.png"/></p>
```
```@raw html
<p align="center"><img width="90%" src="../assets/gui-seq-kspace-load.png"/></p>
```


## Running Simulation
(You can also go to [analog steps using Scripts](programming-workflow.md#Running-Simulation))

Once the inputs are loaded in RAM, it is possible to perform the simulation to get the **Raw Signal**.

### Simulation Parameters

To visualize the default simulation parameters, click on the `Simulate!` dropdown and then press the `View Options` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-sim-params-view.png"/></p>
```

### Visualization of the Raw Signal

Press the `Simulate!` button to perform the simulation (this may take a while). Automatically the generated **Raw Signal** should be displayed or you can click on the `Raw Data` dropdown and then press the `View Raw Data` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-rawsignal-view.png"/></p>
```

## Reconstructing Image using MRIReco
(You can also go to [analog steps using Scripts](programming-workflow.md#Reconstructing-Image-using-MRIReco))

Once the **Raw Signal** is loaded in RAM, it is possible to reconstruct the image.

### Reconstruction Parameters

To visualize the default reconstruction parameters, click on the `Reconstruct!` dropdown and then press the `View Options` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-recon-params-view.png"/></p>
```

### Visualization of the Image

Press the `Reconstruct!` button to perform the reconstruction (this may take a while).  Automatically the generated **Image** should be displayed or you can click on the he `Reconstruct!` dropdown and then press the `|Image|` button:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-image-view.png"/></p>
```

## Exporting Results to .mat File
(You can also go to [analog steps using Scripts](programming-workflow.md#Exporting-Results-to-.mat-File))

The user interface has the option to save the results in **.mat** format. Simply click on the `Export to .mat` and you have the alternatives to get data independently or you can press the `All` button to have all the results given by the simulator:
```@raw html
<p align="center"><img width="90%" src="../assets/gui-export-to-mat.png"/></p>
```

So far, and due to limitations of the user interface dependencies, the **.mat** files are saved in the temporal directory of your computer OS, which can be found by typing the `tempdir()` command in the **Julia REPL**.

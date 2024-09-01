# Getting Started

## Installing Julia

**KomaMRI** was written in Julia, so the first thing you should do is to install it! The latest version of Julia can be downloaded at the [Julia Downloads](https://julialang.org/downloads/) page. It is advisable you add julia to the PATH, which can be done during the installation process.

## Installing KomaMRI

Once Julia is installed, open the Julia REPL, and add the **KomaMRI** package by typing the following commands:

* Press the `]` key and then press `enter` to bring up **Julia**'s package manager.
* Type `add KomaMRI` and then press `enter` in the package manager session.

This process should take about 5 minutes in a fresh Julia installation. Here is how it looks in the **Julia REPL**:

```julia-repl
julia> ]

(@v1.9) pkg> add KomaMRI
```
Then press `Ctrl+C` or `backspace` to return to the `julia>` prompt.


---
## My First MRI Simulation

For our first simulation we will use **KomaMRI**'s graphical user interface (GUI). For this, you will first need to load **KomaMRI** by typing `using KomaMRI`, and then launch the GUI with the [`KomaUI`](@ref) function. Note that if you want to run simulations on the GPU (for example, using CUDA), then `using CUDA` is also necessary (see [GPU Parallelization](../explanation/4-gpu-explanation.md)).  

```julia-repl
julia> using KomaMRI, CUDA

julia> KomaUI()
```
The first time you use this command it may take more time than usual, but a window with the Koma GUI will pop up:

![](../assets/ui-mainpage.png)

The user interface has some basic definitions for the scanner, phantom, and sequence already preloaded. So you can immediately interact with the simulation and reconstruction processes, and then visualize the results.

As a simple demonstration, press the `Simulate!` button and wait until the simulation is ready. Now you have acquired the `Raw Signal` and you should see the following:

![](../assets/ui-view-raw-data.png)

Then, press the `Reconstruct!` button and wait until the reconstruction ends. Now you have reconstructed an `Image` from the `Raw Signal` and you should see the following in the GUI:

![](../assets/ui-view-abs-image.png)

Congratulations, you successfully simulated an MRI acquisition! 🎊

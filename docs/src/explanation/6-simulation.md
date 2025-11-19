# Simulation

## General Overview

**KomaMRI** simulates the magnetization of each spin of a **Phantom** for variable magnetic fields given by a **Sequence**. It is assumed that a single spin is independent of the state of the other spins in the system (a key feature that enables parallelization). Furthermore, there are defined two regimes in the **Sequence**: excitation and precession. During the latter, the excitation fields are nulled and are useful for simplifying some physical equations.

The are more internal considerations in the **KomaMRI** implementation. The **Figure 1** summarizes the functions called to perform the simulation.
```@raw html
<center><img width="100%" src="../../assets/koma-solution.svg"></center>
```
**Figure 1**: The sequence `seq` is discretized after calculating the required time points in the wrapper function [simulate](@ref). The time points are then divided into `Nblocks` to reduce the amount of memory used. The phantom `obj` is divided into `Nthreads`, and **KomaMRI** will use either `run_spin_excitation!` or `run_spin_precession!` depending on the regime. If an [`ADC`](@ref KomaMRIBase.ADC) object is present, the simulator will add the signal contributions of each thread to construct the acquired signal `sig[t]`. All the parameters: `Nthreads`, `Nblocks`, `Δt_rf`, and `Δt`, are passed through a dictionary called `sim_params` as an optional parameter of the simulate function.

From the programming perspective, it is needed to call the [`simulate`](@ref) function with the `sim_params` dictionary keyword argument. A user can change the values of the following keys:

| Parameter | Description |
|:---|:---|
|`"return_type"` | defines the output of the [`simulate`](@ref) function. Possible values are `"raw"`, `"mat"`, and `"state"`, corresponding to outputting a **MRIReco** `RawAcquisitionData`, the signal values, and the last magnetization state of the simulation, respectively. |
| `"sim_method"` | defines the type of simulation. The default value is `Bloch()`, but you can alternatively use the `BlochDict()` simulation method. Moreover, you have the flexibility to create your own methods without altering the **KomaMRI** source code; for further details, refer to the [Simulation Method Extensibility section](#Simulation-Method-Extensibility). |
| `"Δt"` | raster time for gradients. |
| `"Δt_rf"` | raster time for RFs. |
| `"precision"` | defines the floating-point simulation precision. You can choose between `"f32"` and `"f64"` to use `Float32` and `Float64` primitive types, respectively. It's important to note that, especially for GPU operations, using `"f32"` is generally much faster. |
| `"Nblocks"` | divides the simulation into a specified number of time blocks. This parameter is designed to conserve RAM resources, as **KomaMRI** computes a series of simulations consecutively, each with the specified number of blocks determined by the value of `"Nblocks"`. |
| `"Nthreads"` | divides the **Phantom** into a specified number of threads. Because spins are modeled independently of each other, **KomaMRI** can solve simulations in parallel threads, speeding up the execution time. |
| `"gpu"` | is a boolean that determines whether to use GPU or CPU hardware resources, as long as they are available on the host computer. |
| `"gpu_device"` | sets the index ID of the available GPU in the host computer. |

For instance, if you want to perform a simulation on the CPU with float64 precision using the `BlochDict()` method (assuming you have already defined `obj`, `seq` and `sys`), you can do so like this:
```julia
# Set non-default simulation parameters and run simulation
sim_params = KomaMRICore.default_sim_params() 
sim_params["gpu"] = false
sim_params["precision"] = "f64"
sim_params["sim_method"] = BlochDict()
raw = simulate(obj, seq, sys; sim_params)
```

Additionally, the user must be aware of the functions `run_spin_excitation!` and `run_spin_precession!` which defines the algorithm for excitation and precession regimes respectively and can be changed by the user without modifying the source code (more details at [Simulation Method Extensibility](#Simulation-Method-Extensibility)).

Previous simulation, the **Sequence** is discretized to consider specific time points which are critical for simulation. The user can control the time between intermediate gradient samples with the parameter `Δt`. Similarly, the parameter `Δt_rf` manages the time between RF samples, and can be relatively large for 2D imaging where the slice profile is less relevant.

### Computation Efficiency

To reduce the memory usage of our simulator, we subdivided time into `Nblocks`. **KomaMRI** classifies each block in either the excitation regime or the precession regime before the simulation.

We increased the simulation speed by separating the calculations into `Nthreads` and then performing the GPU parallel operations with **CUDA.jl** . This separation is possible as all magnetization vectors are independent of one another.

### Simulation Method Extensibility

In **Julia**, functions use different methods based on the input types via multiple dispatch. We used this to specialize the simulation functions for a given `sim_method <:SimulationMethod` specified in `sim_params`. For a given simulation method, the function `initialize_spin_state` outputs a variable `Xt <: SpinStateRepresentation` that is passed through the simulation (**Figure 1**). For the default simulation method `Bloch`, the spin state is of type `Mag`, but can be extended to a custom representation, like for example EPGs44 or others. Then, the functions `run_spin_excitation!` and `run_spin_precession!` can be described externally for custom types `sim_method` and `Xt`, extending **Koma**’s functionalities without the need of modifying the source code and taking advantage of all of **Koma**’s features.


## Bloch Simulation Method

This is the default simulation method used by **KomaMRI**, however it can always be specified by setting the `sim_method = Bloch()` entry of the `sim_params` dictionary. In the following subsection, we will explain the physical and mathematical background and some considerations and assumptions that enables to speed up the simulation.

### Physical and Mathematical Background

The **Bloch** method of **KomaMRI** simulates the magnetization of each spin by solving the Bloch equations in the rotating frame:
```math
\begin{align} \tag{1}

\frac{\mathrm{d} \boldsymbol{M}}{\mathrm{d} t} =
  \gamma \boldsymbol{M} \times \boldsymbol{B}
- \frac{M_x \hat{x} + M_y \hat{y}}{T_2}
- \frac{M_z \hat{x} + M_0 \hat{y}}{T_1} \:,

\end{align}
```

with ``\gamma`` the gyromagnetic ratio, ``\boldsymbol{M} = [M_x,\: M_y,\: M_z]^T`` the magnetization vector, and
```math
\boldsymbol{B} = [B_{1,x}(t),\: B_{1,y}(t),\: \boldsymbol{G}(t) \cdot \boldsymbol{x} + \Delta \omega(t)]^T
```

the effective magnetic field. ``M_0`` is the proton density, ``T_1`` and ``T_2`` are the relaxation times, and ``\Delta \omega`` is the off-resonance, for each position.

The **Bloch Simulation Method** also uses the technique of **operator splitting** to simplify the solution of Equation `(1)`. This reflects mathematically the intuition of separating the Bloch equations in a rotation operator described by
```math
\begin{align} \tag{2}

\frac{\mathrm{d}}{\mathrm{d}t} \boldsymbol{M} =
\begin{bmatrix}
 0          &  \gamma B_z & -\gamma B_y \\
-\gamma B_z &  0          &  \gamma B_x \\
 \gamma B_y & -\gamma B_x &  0
\end{bmatrix}
\boldsymbol{M} \:,

\end{align}
```

and a relaxation operator described by
```math
\begin{align} \tag{3}

\frac{\mathrm{d}}{\mathrm{d}t} \boldsymbol{M} =
\begin{bmatrix}
-\tfrac{1}{T_2} & 0 & 0 \\
0 & -\tfrac{1}{T_2} & 0 \\
0 & 0 & -\tfrac{1}{T_1}
\end{bmatrix}
\boldsymbol{M}
+
\begin{bmatrix}
0 \\
0 \\
\tfrac{M_0}{T_1}
\end{bmatrix} \:.

\end{align}
```

The evolution of the magnetization can then be described as a two-step process for each time step ``\Delta t`` (**Figure 2**).
```@raw html
<p align="center">
<figure>
  <img width="60%" src="../../assets/block-equation-intuition.svg">
  <figcaption><b>Figure 2</b>: Solution of the Bloch equations for one time step can be described by (2) a rotation and (3) a relaxation step.</figcaption>
</figure>
</p>
```

### Bloch() Method Example

We will consider an RF pulse that excites a phantom with 3 spins, and then we acquire the data:
```@raw html
<details><summary>View code</summary>
```
```julia
# Import modules
using KomaMRI

# Define sequence
ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF, durRF)

nADC = 8192         # number of acquisition samples
durADC = 4000e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay)

seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
```
```@raw html
</details>
```
```julia-repl
julia> obj = Phantom(x=[-0.5e-3; 0.0; 0.5e-3], T1=[1000e-3; 2000e-3; 500e-3], T2=[500e-3; 1000e-3; 2000e-3]);
```
```julia
julia> plot_seq(seq; slider=false)
```
```@raw html
<object type="text/html" data="../../assets/sim-bloch-seq.html" style="width:100%; height:420px;"></object>
```

The resulting signal from the **Bloch()** method is the sum of magnetizations in the transverse plane (x, y):
```julia
# Configure Bloch() simulation method and run simulation
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = Bloch()
sig = simulate(obj, seq, sys; sim_params)
```
```julia-repl
julia> plot(abs.(sig[:,1,1]))
```
```@raw html
<object type="text/html" data="../../assets/sim-bloch-sig.html" style="width:100%; height:420px;"></object>
```

## BlochDict Simulation Method

This is another simulation method defined in the source code of **KomaMRI**. You can specify it by setting the `sim_method = BlochDict()` entry in the `sim_params` dictionary. Additionally, it offers the option to save the resulting signal in the z-component by using `sim_method = BlochDict(save_Mz=true)`. This method allows you to store the magnetizations of all spins in both the transverse plane (x, y) and the longitudinal axis (z) if specified.

### BlochDict() Method Example

We are going to consider the same setup as in the [Bloch() Method Example](#Bloch()-Method-Example). This includes the same excitation, acquisition, and the same 3-spin phantom:
```@raw html
<details><summary>View code</summary>
```
```julia
# Import modules
using KomaMRI, PlotlyJS

# Define sequence
ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF, durRF)

nADC = 8192         # number of acquisition samples
durADC = 4000e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay)

seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
```
```@raw html
</details>
```
```julia-repl
julia> obj = Phantom(x=[-0.5e-3; 0.0; 0.5e-3], T1=[1000e-3; 2000e-3; 500e-3], T2=[500e-3; 1000e-3; 2000e-3]);
```
```julia
julia> plot_seq(seq; slider=false)
```
```@raw html
<object type="text/html" data="../../assets/sim-bloch-seq.html" style="width:100%; height:420px;"></object>
```

The resulting signal from the **BlochDict()** method comprises the individual magnetizations of all spins in both the transverse plane (x, y) and the longitudinal axis (z):
```julia
# Configure BlochDict() simulation method and run simulation
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict(save_Mz=true)
sig = simulate(obj, seq, sys; sim_params)

# Define the plots for traverse and longitudinal magnetizations
pxy = plot(abs.(sig[:,:,1]));
pz = plot(abs.(sig[:,:,2]));
```
```julia-repl
julia> [pxy pz]
```
```@raw html
<object type="text/html" data="../../assets/sim-blochdict-sig.html" style="width:100%; height:420px;"></object>
```
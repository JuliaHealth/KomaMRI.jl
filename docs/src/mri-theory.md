# Simulation Method

## General Overview

**KomaMRI** simulates the magnetization of each spin of a **Phantom** for variable magnetic fields given by a **Sequence**. It is assumed that a single spin is independent of the state of the other spins in the system (a key feature that enables parallelization). Furthermore, there are defined two regimes in the **Sequence**: excitation and precession. During the latter, the excitation fields are nulled and are useful for simplifying some physical equations.

The are more internal considerations in the **KomaMRI** implementation. The **Figure 1** summarizes the functions called to perform the simulation.
```@raw html
<p align="center">
<figure>
  <img width="100%" src="../assets/koma-solution.svg">
  <figcaption><b>Figure 1</b>: The sequence ${\tt seq }$ is discretized after calculating the required time points in the wrapper function ${\tt simulate}$. The time points are then divided into ${\tt Nblocks}$ to reduce the amount of memory used. The phantom ${\tt obj}$ is divided into ${\tt Nthreads}$, and ${\bf KomaMRI}$ will use either ${\tt run\_spin\_excitation!}$ or ${\tt run\_spin\_precession!}$ depending on the regime. If an ${\tt ADC}$ object is present, the simulator will add the signal contributions of each thread to construct the acquired signal ${\tt sig[t]}$. All the parameters: ${\tt Nthreads}$, ${\tt Nblocks}$, ${\tt Δt_{rf}}$, and ${\tt Δt}$, are passed through a dictionary called ${\tt sim\_params}$ as an optional parameter of the simulate function.
</figure>
</p>
```

From the programming perspective, it is needed to call the function ``{\tt simulate()}`` with the ``{\tt sim\_params}`` dictionary keyword argument. A user at least can change the values of the keys:
* ``{\tt Δt}`` and ``{\tt Δt\_{rf}}``, for simulation time refinements,
* ``{\tt Nblocks}`` and ``{\tt Nthreads}``, for computation efficiency, and
* ``{\tt sim\_method}``, for setting different simulation algorithms.

For more details on how to set these parameters while programming in **Julia**, please refer to the [Simulation Parameters Section](programming-workflow.md#Simulation-Parameters).

Additionally, the user must be aware of the functions ``{\tt run\_spin\_excitation!}`` and ``{\tt run\_spin\_precession!}`` which defines the algorithm for excitation and precession regimes respectively and can be changed by the user without modifying the source code (more details at [Simulation Method Extensibility](#Simulation-Method-Extensibility)).

Previous simulation, the **Sequence** is discretized to consider specific time points which are critical for simulation. The user can control the time between intermediate gradient samples with the parameter ``{\tt Δt}``. Similarly, the parameter ``{\tt Δt\_{rf}}`` manages the time between RF samples, and can be relatively large for 2D imaging where the slice profile is less relevant.

### Computation Efficiency

To reduce the memory usage of our simulator, we subdivided time into ``{\tt Nblocks}``. **KomaMRI** classifies each block in either the excitation regime or the precession regime before the simulation.

We increased the simulation speed by separating the calculations into ``{\tt Nthreads}`` and then performing the GPU parallel operations with **CUDA.jl** . This separation is possible as all magnetization vectors are independent of one another.

### Simulation Method Extensibility

In **Julia**, functions use different methods based on the input types via multiple dispatch. We used this to specialize the simulation functions for a given ``{\tt sim\_method <:SimulationMethod}`` specified in ``{\tt sim\_params}``. For a given simulation method, the function ``{\tt initialize\_spin\_state}`` outputs a variable ``{\tt Xt <: SpinStateRepresentation}`` that is passed through the simulation (**Figure 1**). For the default simulation method **Bloch**, the spin state is of type ``{\tt Mag}``, but can be extended to a custom representation, like for example EPGs44 or others. Then, the functions ``{\tt run\_spin\_excitation!}`` and ``{\tt run\_spin\_precession!}`` can be described externally for custom types ``{\tt sim\_method}`` and ``{\tt Xt}``, extending **Koma**’s functionalities without the need of modifying the source code and taking advantage of all of **Koma**’s features.


## Bloch Simulation Method

This is the default simulation method used by **KomaMRI**, however it can always be specified by setting the ``{\tt sim\_method = Bloch()}`` entry of the ``{\tt sim\_params}`` dictionary. In the following subsection, we will explain the physical and mathematical background and some considerations and assumptions that enables to speed up the simulation.

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
  <img width="60%" src="../assets/block-equation-intuition.svg">
  <figcaption><b>Figure 2</b>: Solution of the Bloch equations for one time step can be described by (2) a rotation and (3) a relaxation step.</figcaption>
</figure>
</p>
```

## BlochDict Simulation Method

This is another simulation method defined in the source code of **KomaMRI**. It can be specified by setting the ``{\tt sim\_method = BlochDict()}`` entry of the ``{\tt sim\_params}`` dictionary.

!!! note
    This section is under construction. More explanation of this simulation method is required.


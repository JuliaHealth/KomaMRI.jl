# Simulation Method

## General Overview

**KomaMRI** simulates the magnetization of each spin of a **Phantom** for variable magnetic fields given by a **Sequence**. It is assumed that a single spin is independent of the state of the other spins in the system (a key feature that enables parallelization). Furthermore, there are defined two regimes in the **Sequence**: excitation and precession. During the latter, the excitation fields are nulled and are useful for simplifying some physical equations.

The are more internal considerations in the **KomaMRI** implementation. The **Figure 1** summarizes the functions called to perform the simulation.
```@raw html
<p align="center">
<figure>
  <img width="100%" src="../assets/koma-solution.png">
  <figcaption><b>Figure 1</b>: The sequence ${\tt seq }$ is discretized after calculating the required time points in the wrapper function ${\tt simulate}$. The time points are then divided into ${\tt Nblocks}$ to reduce the amount of memory used. The phantom ${\tt obj}$ is divided into ${\tt Nthreads}$, and ${\bf KomaMRI}$ will use either ${\tt run\_spin\_excitation!}$ or ${\tt run\_spin\_precession!}$ depending on the regime. If an ${\tt ADC}$ object is present, the simulator will add the signal contributions of each thread to construct the acquired signal ${\tt sig[t]}$. All the parameters: ${\tt Nthreads}$, ${\tt Nblocks}$, ${\tt Δt_{rf}}$, and ${\tt Δt}$, are passed through a dictionary called ${\tt sim\_params}$ as an optional parameter of the simulate function.
</figure>
</p>
```

From the programming perspective, it is needed to call the function ``{\tt simulate}`` with the ``{\tt sim\_params}`` dictionary keyword argument. A user at least can change the values of the keys:
* ``{\tt Δt}`` and ``{\tt Δt\_{rf}}``, for simulation time refinements,
* ``{\tt Nblocks}`` and ``{\tt Nthreads}``, for computation efficiency, and
* ``{\tt sim\_method}``, for setting different simulation algorithms.

Additionally, the user must be aware of the functions ``{\tt run\_spin\_excitation!}`` and ``{\tt run\_spin\_precession!}`` which defines the algorithm for excitation and precession regimes respectively and can be changed by the user without modifying the source code (more detail at [Simulation Method Extensibility](#simulation-method-extensibility)).

Previous simulation, the **Sequence** is discretized to consider specific time points which are critical for simulation. The user can control the time between intermediate gradient samples with the parameter ``{\tt Δt}``. Similarly, the parameter ``{\tt Δt\_{rf}}`` manages the time between RF samples, and can be relatively large for 2D imaging where the slice profile is less relevant.

### Computation Efficiency

To reduce the memory usage of our simulator, we subdivided time into ``{\tt Nblocks}``. **KomaMRI** classifies each block in either the excitation regime or the precession regime before the simulation.

We increased the simulation speed by separating the calculations into ``{\tt Nthreads}`` and then performing the GPU parallel operations with **CUDA.jl** . This separation is possible as all magnetization vectors are independent of one another.

### Simulation Method Extensibility

In **Julia**, functions use different methods based on the input types via multiple dispatch. We used this to specialize the simulation functions for a given ``{\tt sim\_method <:SimulationMethod}`` specified in ``{\tt sim\_params}``. For a given simulation method, the function ``{\tt initialize\_spin\_state}`` outputs a variable ``{\tt Xt <: SpinStateRepresentation}`` that is passed through the simulation (**Figure 1**). For the default simulation method **Bloch**, the spin state is of type ``{\tt Mag}``, but can be extended to a custom representation, like for example EPGs44 or others. Then, the functions ``{\tt run\_spin\_excitation!}`` and ``{\tt run\_spin\_precession!}`` can be described externally for custom types ``{\tt sim\_method}`` and ``{\tt Xt}``, extending **Koma**’s functionalities without the need of modifying the source code and taking advantage of all of **Koma**’s features.


## Bloch Simulation Method

This is the default simulation method used by **KomaMRI**, however it can always be specified by setting ``{\tt sim\_method = Bloch()}``. In the following subsection, we will explain the physical and mathematical background and some considerations and assumptions that enables to speed up the simulation.

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

Recall that **KomaMRI** separates the excitation and precession regimes. In the precession regime, the operator splitting method gives an exact solution since the fields ``B_x = B_y = 0`` in Equation `(2)`, whereas during the excitation regime the method has ``O({\Delta t}^3)`` convergence.

From this point forward, we will drop the vectorial notation for ``\boldsymbol{M}`` and ``\boldsymbol{B}_1``, and we will use ``M_{xy} = M_x + i M_y`` and ``B_1 = B_{1,x} + i B_{1,y}`` to describe the simplifications made in each regime.

The rotations during the excitation regime are stored in their spin-domain or SU(2) representation
```math
\bold{Q} =
\begin{bmatrix}
\alpha &-\beta^* \\
\beta  &-\alpha^*
\end{bmatrix}\:, \quad\quad
\text{with}\:
|\alpha|^2 + |\beta|^2 = 1 \:,
```

characterized by the Cayley-Klein complex parameters or Spinors for short ``(\alpha,\:\beta)``. Spinors can represent any 3D
rotation as
```math
\alpha = \cos \left( \tfrac{\varphi}{2} \right)  - i \: n_z \sin \left( \tfrac{\varphi}{2} \right) \\
\beta = -i n_{xy} \sin \left( \tfrac{\varphi}{2} \right) \:.
```

To solve Equation `(2)` the parameters for the Spinors are ``n_{xy} = \tfrac{B_1}{\lVert \boldsymbol{B} \rVert}``, ``n_z = \tfrac{B_z}{\lVert \boldsymbol{B} \rVert}``, and
```math
\begin{align} \tag{4}

\varphi = - \gamma \lVert \boldsymbol{B} \rVert \Delta t \:.

\end{align}
```

Then, the application of a Spinor rotation to a magnetization element is described by the operation
```math
\begin{align} \tag{5}

\begin{bmatrix}
M_{xy}^+ \\
M_z^+
\end{bmatrix} = 
\begin{bmatrix}
2{\alpha}^* \beta M_z + {\alpha^*}^2 M_{xy} - \beta^2 M_{xy}^* \\
(|\alpha|^2 - |\beta|^2)M_z - 2\Re\left( \alpha \beta M_{xy}^* \right)
\end{bmatrix}\:.

\end{align}
```

For the precession regime, all the rotations are with respect to ``z``, and therefore they can be described with a complex exponential applied to the transverse magnetization
```math
\begin{align} \tag{6}

M_{xy}^+ = M_{xy} e^{i\varphi} \:,

\end{align}
```
where ``\varphi`` is defined in Equation `(4)`.

Finally, to solve the relaxation step described in Equation `(3)` the magnetization is updated by
```math
\begin{bmatrix}
M_{xy}^+ \\
M_z^+
\end{bmatrix} =
\begin{bmatrix}
M_{xy} e^{-\tfrac{\Delta t}{T_2}} \\
M_z e^{-\tfrac{\Delta t}{T_1}} + M_0\left(1-e^{-\tfrac{\Delta t}{T_1}}\right)
\end{bmatrix} \:.
```

For precession blocks, we can improve the accuracy of the simulations by using the integral representation of Equation `(6)`, obtained by applying the limit as ``\Delta t \rightarrow 0`` of iterated applications of Equation `(6)`, giving a phase of
```math
\varphi = - \gamma \int_{t_i}^{t_{i+1}} \boldsymbol{G}(\tau) \cdot \boldsymbol{x}(\tau)  \mathrm{d}\tau - \int_{t_i}^{t_{i+1}} \Delta \omega(\tau)  \mathrm{d}\tau \:.
```

Assuming that during the ``i``-th simulation block (``t \in [t_i,\:t_{i+1}]``) the gradients ``\boldsymbol{G}(t)`` are piece-wise linear functions, and ``\boldsymbol{x}(t)`` and ``\Delta \omega (t)`` are approximately constant, then, if we use the trapezoidal rule to obtain the value of this integral, we will obtain an exact result by sampling just the vertices of ``\boldsymbol{G}(t)``, greatly reducing the number of points required by the simulation. We will only need intermediate points in the case of motion and for recording the sampling points as required by the Analog to Digital Converter (ADC). 

We can do something similar with ``B_1(t)`` in the excitation regime. If we assume ``B_1(t)`` is a piece-wise constant function (or concatenation of hard pulses), then Equation `(5)` will give an exact solution to Equation `(2)`. 


## BlochDict Simulation Method

This is another simulation method defined in the source code of **KomaMRI**. It can be specified by setting ``{\tt sim\_method = BlochDict()}``. It performs the same algorithm as **Bloch** method, however it allows you to save the projection of the magnetization in the **z** component.

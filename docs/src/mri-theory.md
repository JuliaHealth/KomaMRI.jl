# Koma's Simulation Method

## Physical and Mathematical Background

Koma simulates the magnetization of each spin by solving the Bloch equations in the rotating frame:
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

The solution of Equation `(1)` for a single spin is independent of the state of the other spins in the system, a key feature that enables parallelization (look at [GPU/CPU Parallelization](#GPU/CPU-Parallelization)).

Our simulator also uses the method of **operator splitting** to simplify the solution of Equation `(1)`. This reflects mathematically the intuition of separating the Bloch equations in a rotation operator described by
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

The evolution of the magnetization can then be described as a two-step process for each time step ``\Delta t`` (**Figure 1**).
```@raw html
<p align="center">
<figure>
  <img width="60%" src="../assets/block-equation-intuition.svg">
  <figcaption><b>Figure 1</b>: Solution of the Bloch equations for one time step can be described by (2) a rotation and (3) a relaxation step.</figcaption>
</figure>
</p>
```

Furthermore, we define two regimes in the pulse sequence: excitation and precession. During the latter, the excitation fields are nulled: ``B_x = B_y = 0`` in Equation `(2)`. In the precession regime, the operator splitting method gives an exact solution, whereas during the excitation regime the method has ``O({\Delta t}^3)`` convergence.

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

## Simulation Blocks, Regime Switching, and Sequence-Aware Time Stepping

To reduce the memory usage of our simulator, we subdivided time into **Nblocks** (**Figure 2**). KomaMRI classifies each block in either the excitation regime or the precession regime before the simulation.

For precession blocks, we can improve the accuracy of the simulations by using the integral representation of Equation `(6)`, obtained by applying the limit as ``\Delta t \rightarrow 0`` of iterated applications of Equation `(6)`, giving a phase of
```math
\varphi = - \gamma \int_{t_i}^{t_{i+1}} \boldsymbol{G}(\tau) \cdot \boldsymbol{x}(\tau)  \mathrm{d}\tau - \int_{t_i}^{t_{i+1}} \Delta \omega(\tau)  \mathrm{d}\tau \:.
```

Assuming that during the ``i``-th simulation block (``t \in [t_i,\:t_{i+1}]``) the gradients ``\boldsymbol{G}(t)`` are piece-wise linear functions, and ``\boldsymbol{x}(t)`` and ``\Delta \omega (t)`` are approximately constant, then, if we use the trapezoidal rule to obtain the value of this integral, we will obtain an exact result by sampling just the vertices of ``\boldsymbol{G}(t)``, greatly reducing the number of points required by the simulation. We will only need intermediate points in the case of motion and for recording the sampling points as required by the Analog to Digital Converter (ADC). The user can control the time between intermediate gradient samples with
the parameter **Δt** (**Figure 2**).

We can do something similar with ``B_1(t)`` in the excitation regime. If we assume ``B_1(t)`` is a piece-wise constant function (or concatenation of hard pulses), then Equation `(5)` will give an exact solution to Equation `(2)`. The parameter **Δt_rf** manages the time between RF samples (**Figure 2**), and can be relatively large for 2D imaging where the slice profile is less relevant.

Thus, **KomaMRI** uses the rationale mentioned above to: (1) call different methods based on the regime of each block, while also (2) obtaining a variable time stepping schedule that adapts to the sequence needs. We named the latter sequence-aware time stepping (**Figure 2**).




## GPU/CPU Parallelization

We further increase the simulation speed by separating the Bloch calculations into **Nthreads** and then performing the GPU operations with CUDA.jl (**Figure 2**). This separation is possible as all magnetization vectors are independent of one another.

```@raw html
<p align="center">
<figure>
  <img width="100%" src="../assets/koma-solution.svg">
  <figcaption><b>Figure 2</b>: This is a summary of the functions called to perform the simulation. The sequence <b>seq</b> is discretized after calculating the required time points in the wrapper function <b>simulate</b>. The time points are then divided into <b>Nblocks</b> to reduce the amount of memory used. The phantom <b>obj</b> is divided into <b>Nthreads</b>, and <b>KomaMRI</b> will use either <b>run_spin_excitation</b> or <b>run_spin_precession</b> depending on the regime. If an ADC object is present, the simulator will add the signal contributions of each thread to construct the acquired signal <b>S[t]</b>. All the parameters: <b>Nthreads</b>, <b>Nblocks</b>, <b>Δt_rf</b>, and <b>Δt</b>, are passed through a dictionary called <b>simParams</b> as an optional parameter of the <b>simulate</b> function.
</figure>
</p>
```

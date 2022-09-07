# MRI Theory

This section is meant to be a general overview or summary of the main MRI concepts and insights. It is a good starting point to show up the most relevant components involved and how they are related for **raw signal** acquisition and image reconstruction. The idea is to have a fresh and clear understanding of what is happening behind the scenes when using the **KomaMRI.jl** package. Some light background in Differential Equations, Signal Processing and Fourier Theory is advisable to follow along.   

## Raw Signal Generation 

In order to generate an **image** from a **phantom** object with a **scanner** system and **sequence** signals, its necessary to acquire a **raw signal** ``s(t)``. This signal can be though as the sum of all spin magnetizations ``M`` of the object:

```math
s(t) =
\int_{x} \int_{y} \int_{z}
\underbrace{
M(x, y, z, t)
}_{\approx \: \alpha \: image(x,y,z)} 
\: \mathrm{d}z \: \mathrm{d}y \: \mathrm{d}x
```

Note that the magnitude of the magnetization is kind of proportional to the **image**. In real life it's not possible to get directly all the spin magnetizations, however it's possible to obtain the sum of all of them in the **raw signal**. To avoid losing image information in the sum operation, every spin magnetization resonates with different Larmor frequencies which values depend on the position ``(x,y,z)`` (i.e. modulated or encoded with the basis of the spatial frequency domain). Thus the **raw signal** can be thought as:

```math
\begin{align} \tag{1}

s(t) \approx
\int_{x} \int_{y} \int_{z}
\underbrace{
m(x, y, z)
}_{\alpha \: image(x,y,z)}
\underbrace{
e^{-j 2 \pi [k_x(t) x + k_y(t) y + k_z(t) z]}
}_{modulation \: basis} 
\: \mathrm{d}z \: \mathrm{d}y \: \mathrm{d}x

\end{align}
```

where:

```math
\vec{k}(t) =
\begin{pmatrix}
k_x(t) \\
k_y(t) \\
k_z(t)
\end{pmatrix} =
\frac{2\pi}{\gamma}
\begin{pmatrix}
\int_{0}^{t} G_x(\tau) \mathrm{d} \tau\\
\int_{0}^{t} G_y(\tau) \mathrm{d} \tau\\
\int_{0}^{t} G_z(\tau) \mathrm{d} \tau\\
\end{pmatrix}
\;\;\; , \;\;\;
\begin{matrix*}[l]
\gamma: \: gyromagnetic \: ratio \\
G_i(t): \: input \: gradient \: signals
\end{matrix*}
```

In the above expressions, we can see that the frequency of each spin can be manipulated by applying input **gradient** signals (or a gradient field). In practice, this gradient field is applied in the longitudinal axis ``\hat{z}`` (but it is always dependent on the ``(x,y,z)`` position), which makes the spins able to resonate (or precess) at different Larmor frequencies after another input **RF** pulse excite them in the transverse direction. Both inputs, the **gradient** signal and the **RF** pulse, are part of the effective magnetic field ``\vec{B}(t)``:

```math
\vec{B}(t) = 
\begin{pmatrix}
B_{1,x}(t) \\
B_{1,y}(t) \\
G_x(t) x + G_y(t) y + G_z(t) z\\
\end{pmatrix}
\;\;\; , \;\;\;
\begin{matrix*}[l]
B_{1,i}(t): \: input \: RF \: pulse \: (transverse) \\
G_i(t):     \: input \: gradients \: (longitudinal)
\end{matrix*}
```

It's important to highlight that the coil that senses the **raw signal** can only detects magnetization components oriented in the transverse direction. For this reason is necessary to apply the short **RF** signal orthogonally to the longitudinal ``\hat{z}`` axe.

One of the primary concerns, to generate an image is to design proper input signals for the effective magnetic field ``\vec{B}(t)``. In particular, by inspecting equation `(1)`, it's possible to manipulate the **spacial frequencies** ``k_x(t)``, ``k_y(t)`` and ``k_z(t)`` by applying the **gradients** ``G_x(t)``, ``G_y(t)`` and ``G_z(t)``. Thus, we have information of the **raw signal** ``s(t)`` an the basis ``e^{-j 2 \pi [k_x(t) x + k_y(t) y + k_z(t) z]}``. Mathematically speaking, every sample of the **raw signal** is the Fourier transform of the magnetization ``m(x,y,z)`` for a specific point of the ``spacial frequency`` domain:

```math
s(t) = Fourier\{\:m \: \}(k_x(t),\: k_y(t),\: k_z(t))
```

Therefore, to get the magnetization ``m(x,y,z)`` for all the points in the ``spacial`` domain its necessary to solve the inverse problem with enough points to cover the complete ``spacial frequency`` domain, which can be achieved by following a trajectory over time applying different gradient signals (i.e. a trajectory to complete the **k-space**).


## K-Space and Acquisition

Note that the trajectory to cover the **k-space** eventually can have any continuos shape, however it cannot fill the complete space. Furthermore, due to natural hardware restrictions, the continuos trajectory is sampled during the acquisition of the **raw signal** ``s[t]``. Thus, every discrete point of ``s[t]`` represents a discrete point in the **k-space**.

Intuitively, it is desirable to get many points as possible and homogeneously distributed in the **k-space**. In particular, since the theory behind the **raw signal generation** is intimately related with the Fourier Transform, a natural way to cover the **k-space** is by taken a discrete mesh grid of points (trajectories and samples separated by small cubes). In this case, it is possible to apply Fourier theory to define the minimal **k-space** resolution (separation of the samples in the **k-space**) to set space dimensions (Field of Views) and prevent aliasing in the image, and define maximal limits in the **k-space** to set space resolution in the image.

```math
\underbrace{
\Delta k_i
}_{k-space \: resolution}
\longrightarrow
\underbrace{
FOV_i
}_{space \: width \: bounds} 
```
```math
\underbrace{
W_{k_i}
}_{k-space \: width \: bounds}
\longrightarrow
\underbrace{
\Delta i
}_{space \: resolution} 
```
```math
\forall i \in \{x,\: y, \: z\}
```

Even though a mesh grid of discrete points is the natural way of thinking to cover the **k-space**, it is always possible possible to apply more exotic **k-space** trajectories, which could be helpful, for instance, to reduce the complete acquisition time. Keep in mind though, this fact must be regarded when solving the inverse problem for obtaining the image, for example by applying and interpolation function before taking the inverse Fourier Transform.


## Spin Dynamics

It's important to point out that all the magnetization spins are independent from each other, so we could separate the **phantom** object into multiple spins and solve the **Bloch Equations** for every magnetization vector ``\vec{M}`` independently:

```math
\frac{\mathrm{d} \vec{M}}{\mathrm{d} t} =
  \gamma \vec{M} \times \vec{B}
- \frac{M_x \hat{x} + M_y \hat{y}}{T_2}
- \frac{M_z \hat{x} + M_0 \hat{y}}{T_1}
```

or:

```math
\begin{align} \tag{2}

\frac{\mathrm{d}}{\mathrm{d}t} \vec{M} =
\underbrace{
\gamma
\begin{bmatrix}
 0   &  B_z & -B_y \\
-B_z &  0   &  B_x \\
 B_y & -B_x &  0
\end{bmatrix}
\vec{M}
}_\text{rotation} 
-
\underbrace{
\begin{bmatrix}
\tfrac{1}{T_2} & 0 & 0 \\
0 & \tfrac{1}{T_2} & 0 \\
0 & 0 & \tfrac{1}{T_1}
\end{bmatrix}
\vec{M}
}_\text{relaxation} 
+
\underbrace{
\begin{bmatrix}
0 \\
0 \\
\tfrac{M_0}{T_1}
\end{bmatrix}
}_\text{steady-state} 

\end{align}
```

```math
\begin{matrix*}[l]
\gamma: & gyromagnetic \: ratio \\
T_2:    & transverse \: relaxation \: time \: constant \\
T_1:    & longitudinal \: relaxation \: time \: constant
\end{matrix*}
```

```math
\vec{M}(t) =
\begin{pmatrix}
M_x(t) \\
M_y(t) \\
M_z(t)
\end{pmatrix}
\;\;\; , \;\;\;
\vec{B}(t) = 
\begin{pmatrix}
B_x(t) \\
B_y(t) \\
B_z(t)
\end{pmatrix} =
\begin{pmatrix}
B_{1,x}(t) \\
B_{1,y}(t) \\
G_x(t) x + G_y(t) y + G_z(t) z\\
\end{pmatrix}
```

```math
\begin{matrix*}[l]
B_{1,i}(t): & input \: RF \: pulse \: (transverse) \\
G_i(t):     & input \: gradients \: (longitudinal)
\end{matrix*}
```

Note that equation `(2)` can be separated into three parts:

* Rotation: governed by the inputs RF pulse and gradient signals. It gives an initial excitation and the oscillatory behavior for different Larmor frequencies, respectively. 
* Relaxation: gives the decay behavior (the magnetization envelope) after the excitation of the spins.
* Steady-State: spins points towards the longitudinal direction after a while.

Thus, to get the **raw signal** ``s(t)`` it's necessary to solve the **Bloch equations** (equation `(2)`) for every spin of the **phantom** object, then sum up the contributions of all of them and finally consider just the components of the transverse plane:

```math
s(t) = s_{xy}(t)
\;\;\; , \;\;\;
s_{xy}(t) = s_{x}(t) + j s_{y}(t)
```

```math
\begin{pmatrix}
s_{x}(t) \\
s_{y}(t) \\
s_{z}(t)
\end{pmatrix} =
\int_{x} \int_{y} \int_{z}
\vec{M}(x, y, z, t)
\: \mathrm{d}z \: \mathrm{d}y \: \mathrm{d}x
```
abstract type BlochMagnus <: SimulationMethod end

@doc raw"""
    BlochMagnusConst1()

Bloch-Magnus simulation with a piecewise-constant effective angular frequency
``\boldsymbol{\omega}(\tau)`` [rad/s] and the first Magnus term.

Expected smooth global error: ``\mathcal{O}(\Delta t)``.

Alias: `BlochMagnus1()`.

```math
\boldsymbol{\theta}_1 = \boldsymbol{\omega}_0 \Delta t
```
"""
struct BlochMagnusConst1 <: BlochMagnus end

@doc raw"""
    BlochMagnusLin2()

Bloch-Magnus simulation with piecewise-linear ``\boldsymbol{\omega}(\tau)``
[rad/s] samples ``\boldsymbol{\omega}_0`` and ``\boldsymbol{\omega}_1``.

Expected smooth global error: ``\mathcal{O}(\Delta t^2)``.

```math
\boldsymbol{\theta}_1 =
\frac{\Delta t}{2}(\boldsymbol{\omega}_0 + \boldsymbol{\omega}_1)
```
"""
struct BlochMagnusLin2 <: BlochMagnus end

@doc raw"""
    BlochMagnusMid2()

Bloch-Magnus simulation with one midpoint sample
``\boldsymbol{\omega}_{1/2}`` [rad/s].

Expected smooth global error: ``\mathcal{O}(\Delta t^2)``.

Alias: `BlochMagnus2()`.

```math
\boldsymbol{\theta}_1 =
\Delta t\,\boldsymbol{\omega}_{1/2}
```
"""
struct BlochMagnusMid2 <: BlochMagnus end

@doc raw"""
    BlochMagnusLinComm2()

Bloch-Magnus simulation with piecewise-linear ``\boldsymbol{\omega}(\tau)``
[rad/s] samples ``\boldsymbol{\omega}_0`` and ``\boldsymbol{\omega}_1``,
including the endpoint commutator correction.

Expected smooth global error: ``\mathcal{O}(\Delta t^2)``.

```math
\boldsymbol{\theta}_1 =
\frac{\Delta t}{2}(\boldsymbol{\omega}_0 + \boldsymbol{\omega}_1)
```

```math
\boldsymbol{\theta}_2 =
-\frac{\Delta t^2}{12}
(\boldsymbol{\omega}_0 \times \boldsymbol{\omega}_1)
```
"""
struct BlochMagnusLinComm2 <: BlochMagnus end

@doc raw"""
    BlochMagnusQuad2()

Bloch-Magnus simulation with piecewise-quadratic ``\boldsymbol{\omega}(\tau)``
[rad/s] samples ``\boldsymbol{\omega}_0``, ``\boldsymbol{\omega}_{1/2}``,
and ``\boldsymbol{\omega}_1``.

Expected smooth global error: ``\mathcal{O}(\Delta t^2)``.

```math
\boldsymbol{\theta}_1 =
\frac{\Delta t}{6}
(\boldsymbol{\omega}_0 + 4\boldsymbol{\omega}_{1/2} + \boldsymbol{\omega}_1)
```
"""
struct BlochMagnusQuad2 <: BlochMagnus end

@doc raw"""
    BlochMagnusQuad4()

Bloch-Magnus simulation with piecewise-quadratic ``\boldsymbol{\omega}(\tau)``
[rad/s] samples ``\boldsymbol{\omega}_0``, ``\boldsymbol{\omega}_{1/2}``,
and ``\boldsymbol{\omega}_1``.

Expected smooth global error: ``\mathcal{O}(\Delta t^4)``.

```math
\boldsymbol{\theta}_1 =
\frac{\Delta t}{6}
(\boldsymbol{\omega}_0 + 4\boldsymbol{\omega}_{1/2} + \boldsymbol{\omega}_1)
```

Let

```math
\boldsymbol{\delta} =
\boldsymbol{\omega}_{1/2} -
\frac{\boldsymbol{\omega}_0 + \boldsymbol{\omega}_1}{2}.
```

Then

```math
\boldsymbol{\theta}_2 =
-\frac{\Delta t^2}{12}
(\boldsymbol{\omega}_0 \times \boldsymbol{\omega}_1)
+
\frac{\Delta t^2}{15}
\left((\boldsymbol{\omega}_1 - \boldsymbol{\omega}_0) \times \boldsymbol{\delta}\right)
```
"""
struct BlochMagnusQuad4 <: BlochMagnus end

@doc raw"""
    BlochMagnusGL2()

Bloch-Magnus simulation with two Gauss-Legendre samples
``\boldsymbol{\omega}_-`` and ``\boldsymbol{\omega}_+`` [rad/s].

Expected smooth global error: ``\mathcal{O}(\Delta t^2)``.

```math
\tau = \frac{t - t_n}{\Delta t}, \qquad
\boldsymbol{\omega}_- =
\boldsymbol{\omega}\left(\frac{1}{2} - \frac{\sqrt{3}}{6}\right),
\qquad
\boldsymbol{\omega}_+ =
\boldsymbol{\omega}\left(\frac{1}{2} + \frac{\sqrt{3}}{6}\right).
```

```math
\boldsymbol{\theta}_1 =
\frac{\Delta t}{2}(\boldsymbol{\omega}_- + \boldsymbol{\omega}_+)
```
"""
struct BlochMagnusGL2 <: BlochMagnus end

@doc raw"""
    BlochMagnusGL4()

Bloch-Magnus simulation with two Gauss-Legendre samples
``\boldsymbol{\omega}_-`` and ``\boldsymbol{\omega}_+`` [rad/s].

Expected smooth global error: ``\mathcal{O}(\Delta t^4)``.

Alias: `BlochMagnus4()`.

```math
\tau = \frac{t - t_n}{\Delta t}, \qquad
\boldsymbol{\omega}_- =
\boldsymbol{\omega}\left(\frac{1}{2} - \frac{\sqrt{3}}{6}\right),
\qquad
\boldsymbol{\omega}_+ =
\boldsymbol{\omega}\left(\frac{1}{2} + \frac{\sqrt{3}}{6}\right).
```

```math
\boldsymbol{\theta}_1 =
\frac{\Delta t}{2}(\boldsymbol{\omega}_- + \boldsymbol{\omega}_+)
```

```math
\boldsymbol{\theta}_2 =
-\frac{\sqrt{3}\Delta t^2}{12}
(\boldsymbol{\omega}_- \times \boldsymbol{\omega}_+)
```
"""
struct BlochMagnusGL4 <: BlochMagnus end

@doc raw"""
    BlochMagnusBGL4()

Bloch-Magnus simulation using the fourth-order Blanes Gauss-Legendre method.

Expected smooth global error: ``\mathcal{O}(\Delta t^4)``.

It samples ``\boldsymbol{\omega}(\tau)`` [rad/s] at
``\boldsymbol{\omega}_-``, ``\boldsymbol{\omega}_{1/2}``, and
``\boldsymbol{\omega}_+``.

Let

```math
\tau = \frac{t - t_n}{\Delta t}, \qquad
c = \sqrt{\frac{3}{20}}, \qquad
\boldsymbol{\omega}_- =
\boldsymbol{\omega}\left(\frac{1}{2} - c\right),
```

```math
\boldsymbol{\omega}_{1/2} =
\boldsymbol{\omega}\left(\frac{1}{2}\right), \qquad
\boldsymbol{\omega}_+ =
\boldsymbol{\omega}\left(\frac{1}{2} + c\right).
```

Define

```math
\boldsymbol{i}_0 =
\frac{1}{18}
(5\boldsymbol{\omega}_- + 8\boldsymbol{\omega}_{1/2} + 5\boldsymbol{\omega}_+),
\qquad
\boldsymbol{i}_1 =
\frac{\sqrt{15}}{36}(\boldsymbol{\omega}_+ - \boldsymbol{\omega}_-).
```

Then

```math
\boldsymbol{\theta}_1 = \Delta t\,\boldsymbol{i}_0
```

and

```math
\boldsymbol{\theta}_2 =
-\Delta t^2(\boldsymbol{i}_0 \times \boldsymbol{i}_1).
```

# References
- Blanes, Casas, Oteo, and Ros,
  [The Magnus expansion and some of its applications](https://doi.org/10.1016/j.physrep.2008.11.001),
  Physics Reports 470 (2009).
- Stephanie Gonzalez, *Solving the Bloch equation with the Magnus expansion*
  (master's thesis).
"""
struct BlochMagnusBGL4 <: BlochMagnus end

@doc raw"""
    BlochMagnusBGL6()

Bloch-Magnus simulation using the sixth-order Blanes Gauss-Legendre method.

Alias: `BlochMagnus6()`.

Expected smooth global error: ``\mathcal{O}(\Delta t^6)``.

It samples ``\boldsymbol{\omega}(\tau)`` [rad/s] at
``\boldsymbol{\omega}_-``, ``\boldsymbol{\omega}_{1/2}``, and
``\boldsymbol{\omega}_+``.

Let

```math
\tau = \frac{t - t_n}{\Delta t}, \qquad
c = \sqrt{\frac{3}{20}}, \qquad
\boldsymbol{\omega}_- =
\boldsymbol{\omega}\left(\frac{1}{2} - c\right),
```

```math
\boldsymbol{\omega}_{1/2} =
\boldsymbol{\omega}\left(\frac{1}{2}\right), \qquad
\boldsymbol{\omega}_+ =
\boldsymbol{\omega}\left(\frac{1}{2} + c\right).
```

Define

```math
\boldsymbol{i}_0 =
\frac{1}{18}
(5\boldsymbol{\omega}_- + 8\boldsymbol{\omega}_{1/2} + 5\boldsymbol{\omega}_+),
\qquad
\boldsymbol{i}_1 =
\frac{\sqrt{15}}{36}(\boldsymbol{\omega}_+ - \boldsymbol{\omega}_-),
```

```math
\boldsymbol{i}_2 =
\frac{1}{24}(\boldsymbol{\omega}_- + \boldsymbol{\omega}_+).
```

Then, with
``\boldsymbol{j} = \frac{3}{2}\boldsymbol{i}_0 - 6\boldsymbol{i}_2``
and
``\operatorname{box}(\boldsymbol{u}, \boldsymbol{v}) =
\boldsymbol{u} \times (\boldsymbol{u} \times \boldsymbol{v}) =
\boldsymbol{u}(\boldsymbol{u} \cdot \boldsymbol{v}) -
\boldsymbol{v}(\boldsymbol{u} \cdot \boldsymbol{u})``,

```math
\boldsymbol{\theta}_1 = \Delta t\,\boldsymbol{i}_0,
```

```math
\boldsymbol{\theta}_2 =
-\Delta t^2(\boldsymbol{j} \times \boldsymbol{i}_1),
```

```math
\boldsymbol{\theta}_{34} =
\Delta t^2\,\operatorname{box}\left(
\boldsymbol{i}_0,
\frac{1}{2}\Delta t\,\boldsymbol{i}_2 - \frac{1}{60}\boldsymbol{\theta}_2
\right)
+
\frac{3}{5}\Delta t^3\,\operatorname{box}(\boldsymbol{i}_1, \boldsymbol{j}),
```

and

```math
\boldsymbol{\theta} =
\boldsymbol{\theta}_1 + \boldsymbol{\theta}_2 + \boldsymbol{\theta}_{34}.
```

# References
- Blanes, Casas, Oteo, and Ros,
  [The Magnus expansion and some of its applications](https://doi.org/10.1016/j.physrep.2008.11.001),
  Physics Reports 470 (2009).
- Stephanie Gonzalez, *Solving the Bloch equation with the Magnus expansion*
  (master's thesis).
"""
struct BlochMagnusBGL6 <: BlochMagnus end

const BlochMagnus1 = BlochMagnusConst1
const BlochMagnus2 = BlochMagnusMid2
const BlochMagnus4 = BlochMagnusGL4
const BlochMagnus6 = BlochMagnusBGL6

export BlochMagnus, BlochMagnusConst1, BlochMagnusLin2, BlochMagnusMid2, BlochMagnusLinComm2, BlochMagnusQuad2, BlochMagnusQuad4, BlochMagnusGL2, BlochMagnusGL4, BlochMagnusBGL4, BlochMagnusBGL6
export BlochMagnus1, BlochMagnus2, BlochMagnus4, BlochMagnus6
const BlochLikeSimMethods = Union{
    Bloch,
    BlochMagnusConst1,
    BlochMagnusLin2,
    BlochMagnusMid2,
    BlochMagnusLinComm2,
    BlochMagnusQuad2,
    BlochMagnusQuad4,
    BlochMagnusGL2,
    BlochMagnusGL4,
    BlochMagnusBGL4,
    BlochMagnusBGL6,
}

include("IntegrationNodes.jl")
include("sampling/IntegrationNodeSampling.jl")

include("rotation_vectors/MagnusConst1.jl")
include("rotation_vectors/MagnusLin2.jl")
include("rotation_vectors/MagnusMid2.jl")
include("rotation_vectors/MagnusLinComm2.jl")
include("rotation_vectors/MagnusQuad2.jl")
include("rotation_vectors/MagnusQuad4.jl")
include("rotation_vectors/MagnusGL2.jl")
include("rotation_vectors/MagnusGL4.jl")
include("rotation_vectors/MagnusBGL4.jl")
include("rotation_vectors/MagnusBGL6.jl")

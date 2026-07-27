# Use Custom Coil Sensitivity Maps

Use `ArbitraryCoilSens` to simulate measured or calculated complex receive-coil
maps. KomaMRI expects the map dimensions in `(x, y, z, coil)` order and the grid
coordinates in metres.

For 2D maps shaped `(Nx, Ny, ncoils)`, add a singleton `z` dimension and tile it
over a short physical range:

```julia
using KomaMRI

Nx, Ny, ncoils = size(maps_2d)
x = range(-FOVx / 2, FOVx / 2; length=Nx)
y = range(-FOVy / 2, FOVy / 2; length=Ny)
z = range(-1e-3, 1e-3; length=3)

maps = repeat(
    reshape(maps_2d, Nx, Ny, 1, ncoils),
    1, 1, length(z), 1,
)
receiver = ArbitraryCoilSens(x, y, z, maps)
sys = Scanner(; receiver)
raw = simulate(obj, seq, sys)
```

The tiling repeats each coil's 2D map only along `z`; it does not duplicate or
combine coil channels. During simulation, KomaMRI linearly interpolates these
samples at the phantom positions. Samples outside the supplied grid have zero
sensitivity.

## Compare direct and SENSE reconstructions

The following measured and simulated two-fold accelerated acquisitions use
ESPIRiT maps estimated from a fully sampled reference acquisition. The direct
reconstruction retains the expected fold-over, while SENSE combines the coil
channels with the complex maps to recover one unfolded image. Hover, pan, or
scroll to inspect each reconstruction interactively.

### Measured acquisition

```@raw html
<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem;">
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/acquisition_direct.html" style="width:100%;height:520px;"></object>
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/acquisition_sense.html" style="width:100%;height:520px;"></object>
</div>
```

### Simulated MRD

```@raw html
<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem;">
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/simulated_mrd_direct.html" style="width:100%;height:520px;"></object>
<object type="text/html" data="../assets/use-custom-coil-sensitivity-maps/simulated_mrd_sense.html" style="width:100%;height:520px;"></object>
</div>
```

Ensure the map grid covers the complete phantom. A masked or undersized map is
zero outside its support and therefore removes signal there during simulation.

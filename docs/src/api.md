# API Documentation

```@contents
Pages = ["api.md"]
Depth = 3
```

## Structs

### `Scanner`
```@docs
Scanner
```

### `Phantom`
```@docs
Phantom
brain_phantom2D
brain_phantom3D
```

### `Sequence`
```@docs
Sequence
```

### `Grad`
```@docs
Grad
Grad(::Function, ::Real, ::Int64)
```

### `RF`
```@docs
RF
```

### `ADC`
```@docs
ADC
```

### `Delay`
```@docs
Delay
```

## Read Data

### `read_seq`
```@docs
KomaMRI.read_seq
```

### `read_phantom_jemris`
```@docs
read_phantom_jemris
```

### `rawSignalToISMRMRD`
```@docs
rawSignalToISMRMRD
```

## Pulse Design

### `PulseDesigner`
```@docs
PulseDesigner
```

### `PulseDesigner.RF_hard`
```@docs
PulseDesigner.RF_hard
```

### `PulseDesigner.EPI`
```@docs
PulseDesigner.EPI
```

### `PulseDesigner.radial_base`
```@docs
PulseDesigner.radial_base
```

## Simulation

### `simulate`
```@docs
simulate
```

### `simulate_slice_profile`
```@docs
simulate_slice_profile
```

## Plots

### `plot_phantom_map`
```@docs
plot_phantom_map
```

### `plot_seq`
```@docs
plot_seq
```

### `plot_kspace`
```@docs
plot_kspace
```

### `plot_signal`
```@docs
plot_signal
```

### `plot_M0`
```@docs
plot_M0
```

### `plot_image`
```@docs
plot_image
```

# KomaMRIBase

```@meta
CurrentModule = KomaMRIBase
```

## `Scanner`-related functions

```@docs
Scanner
```

## `Phantom`-related functions 

```@docs
Phantom
brain_phantom2D
brain_phantom3D
pelvis_phantom2D
heart_phantom
```

### `MotionModel`-related functions

```@docs
get_spin_coords
```

### `SimpleMotion <: MotionModel`

```@docs
SimpleMotion
```

### `SimpleMotion types`

```@docs
Translation
Rotation
HeartBeat
PeriodicTranslation
PeriodicRotation
PeriodicHeartBeat
```

### `ArbitraryMotion <: MotionModel`

```@docs
ArbitraryMotion
```

## `Sequence`-related functions

```@docs
Sequence
dur
get_block_start_times
get_flip_angles
```

### `Grad`

```@docs
Grad
Grad(::Function, ::Real, ::Int64)
```
### `RF`

```@docs
RF
RF(::Function, ::Real, ::Int64)
get_flip_angle
```

### `ADC`

```@docs
ADC
get_adc_sampling_times
get_adc_phase_compensation
```

### `Delay`
    
```@docs
Delay
```

### Rotation matrices

```@docs
rotx
roty
rotz
```

### Moments

```@docs
get_Mk
get_kspace
get_M0
get_M1
get_M2
```

### Event checks

```@docs
is_RF_on
is_GR_on
is_Gx_on
is_Gy_on
is_Gz_on
is_ADC_on
```

### `DiscreteSequence`

```@docs
DiscreteSequence
discretize
get_samples
times
ampls
freqs
```

### Other functions

```@docs
trapz
cumtrapz
kfoldperm
```

## Sequence Building Blocks (SBB)

```@docs
PulseDesigner
PulseDesigner.RF_hard
PulseDesigner.RF_sinc
PulseDesigner.EPI
PulseDesigner.radial_base
PulseDesigner.spiral_base
PulseDesigner.EPI_example
```
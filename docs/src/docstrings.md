# API Documentation

## [Dataflow Graph](@id dataflow-graph)

![](assets/dataflow.svg)

```@contents
Pages = ["docstrings.md"]
Depth = 3
```

## [DataTypes](@id datatypes)

### `Mag`
```@docs
Mag
Mag(::Phantom, ::Symbol)
show(::IO, ::Mag)
getproperty(::Vector{Mag}, ::Symbol)
```

### `Phantom`

Refer to [API Documentation](api.md#Phantom):
* [`KomaMRI.Phantom`](@ref) — Type
* [`KomaMRI.brain_phantom2D`](@ref) — Function
* [`KomaMRI.brain_phantom3D`](@ref) — Function

```@docs
KomaMRI.heart_phantom
```

### `Scanner`

Refer to [API Documentation](api.md#Scanner):
* [`KomaMRI.Scanner`](@ref) — Type

### `Sequence`

Refer to [API Documentation](api.md#Sequence):
* [`KomaMRI.Sequence`](@ref) — Type

```@docs
show(::IO, ::Sequence)
KomaMRI.is_ADC_on
KomaMRI.is_RF_on
KomaMRI.is_GR_on
KomaMRI.is_Gx_on
KomaMRI.is_Gy_on
KomaMRI.is_Gz_on
KomaMRI.is_Delay
KomaMRI.durs(::Sequence)
dur(::Sequence)
KomaMRI.⏢
KomaMRI.get_grads
KomaMRI.get_rfs
KomaMRI.get_flip_angles
KomaMRI.get_ADC_on
KomaMRI.get_kspace
KomaMRI.get_RF_types

KomaMRI.δ2N
```

### `Grad`

Refer to [API Documentation](api.md#Grad):
* [`KomaMRI.Grad`](@ref) — Type
* [`KomaMRI.Grad(::Function, ::Real, ::Int64)`](@ref) — Method

```@docs
rotx
roty
rotz
show(::IO, ::Grad)
getproperty(::Vector{Grad}, ::Symbol)
dur(::Grad)
```

### `RF`

Refer to [API Documentation](api.md#RF):
* [`KomaMRI.RF`](@ref) — Type

```@docs
Spinor
show(::IO,::Spinor)
*(::Spinor, ::Spinor)
Rz
Ry
Rx
KomaMRI.Rg
KomaMRI.Rφ
Q
abs(::Spinor)

show(::IO, ::RF)
getproperty(::Vector{RF}, ::Symbol)
dur(::RF)
KomaMRI.RF_fun
KomaMRI.get_flip_angle
KomaMRI.get_RF_center
```

### `ADC`

Refer to [API Documentation](api.md#ADC):
* [`KomaMRI.ADC`](@ref) — Type

```@docs
getproperty(::Vector{ADC}, ::Symbol)
KomaMRI.get_sample_times
KomaMRI.get_sample_phase_compensation
```

### `Delay`

Refer to [API Documentation](api.md#Delay):
* [`KomaMRI.Delay`](@ref) — Type

```@docs
show(::IO, ::Delay)
+(::Sequence, ::Delay)
```

## [Pulseq.jl](@id pulseq)

Refer to [API Documentation](api.md#read_seq):
* [`KomaMRI.read_seq`](@ref) — Function

### `read_Grad`
```@docs
KomaMRI.read_Grad
```

### `read_RF`
```@docs
KomaMRI.read_RF
```

### `read_ADC`
```@docs
KomaMRI.read_ADC
```

### `get_block`
```@docs
KomaMRI.get_block
```

## [JEMRIS.jl](@id jemris)

### `read_phantom_jemris`
Refer to [API Documentation](api.md#read_phantom_jemris):
* [`KomaMRI.read_phantom_jemris`](@ref) — Function

## [ISMRMRD.jl](@id ismrmrd)

### `rawSignalToISMRMRD`
Refer to [API Documentation](api.md#rawSignalToISMRMRD):
* [`KomaMRI.rawSignalToISMRMRD`](@ref) — Function

## [PulseDesigner.jl](@id pulse-designer)

Refer to [API Documentation](api.md#Pulse-Design):
* [`KomaMRI.PulseDesigner`](@ref) — Function
* [`KomaMRI.PulseDesigner.RF_hard`](@ref) — Function
* [`KomaMRI.PulseDesigner.EPI`](@ref) — Function
* [`KomaMRI.PulseDesigner.radial_base`](@ref) — Function

## [KeyValuesCalculation.jl](@id key-values-calculation)

### `get_theo_A`
```@docs
KomaMRI.get_theo_A
```

### `get_theo_t`
```@docs
KomaMRI.get_theo_t
```

### `get_theo_Gi`
```@docs
KomaMRI.get_theo_Gi
```

## [TrapezoidalIntegration.jl](@id trapezoidal-integration)

### `trapz`
```@docs
KomaMRI.trapz
```

### `cumtrapz`
```@docs
KomaMRI.cumtrapz
```

## [TimeStepCalculation.jl](@id time-step-calculation)

### `points_from_key_times`
```@docs
KomaMRI.points_from_key_times
```

### `get_variable_times`
```@docs
KomaMRI.get_variable_times
```

### `get_uniform_times`
```@docs
KomaMRI.get_uniform_times
```

### `kfoldperm`
```@docs
KomaMRI.kfoldperm
```

### `get_breaks_in_RF_key_points`
```@docs
KomaMRI.get_breaks_in_RF_key_points
```

## [SimulationCore.jl](@id simulation-core)

Refer to [API Documentation](api.md#Simulation):
* [`KomaMRI.simulate`](@ref) — Function
* [`KomaMRI.simulate_slice_profile`](@ref) — Function

### `print_gpus`
```@docs
KomaMRI.print_gpus
```

### `run_spin_precession`
```@docs
KomaMRI.run_spin_precession
```

### `run_spin_precession_parallel`
```@docs
KomaMRI.run_spin_precession_parallel
```

### `run_spin_excitation`
```@docs
KomaMRI.run_spin_excitation
```

### `run_spin_excitation_parallel`
```@docs
KomaMRI.run_spin_excitation_parallel
```

### `run_sim_time_iter`
```@docs
KomaMRI.run_sim_time_iter
```

## [DisplayFunctions.jl](@id display-functions)

Refer to [API Documentation](api.md#Plots):
* [`KomaMRI.plot_seq`](@ref) — Function
* [`KomaMRI.plot_image`](@ref) — Function
* [`plot_kspace`](@ref) — Function
* [`plot_M0`](@ref) — Function
* [`plot_phantom_map`](@ref) — Function
* [`plot_signal`](@ref) — Function

### `theme_chooser`
```@docs
KomaMRI.theme_chooser
```

### `interp_map`
```@docs
KomaMRI.interp_map
```

### `plot_dict`
```@docs
KomaMRI.plot_dict
```

# Home

## [Dataflow Graph](@id dataflow-graph)

![](assets/dataflow.svg)

```@contents
Pages = ["software.md"]
Depth = 3
```

## [datatypes](@id datatypes)

### `Mag`
```@docs
Mag
```

### `Phantom`
```@docs
Phantom
KomaMRI.heart_phantom
brain_phantom2D
brain_phantom3D
```

### `Scanner`
```@docs
Scanner
```

### `Sequence`
```@docs
Sequence
KomaMRI.is_ADC_on
KomaMRI.is_RF_on
KomaMRI.is_GR_on
KomaMRI.is_Gx_on
KomaMRI.is_Gy_on
KomaMRI.is_Gz_on
KomaMRI.is_Delay
KomaMRI.⏢
KomaMRI.get_grads
KomaMRI.get_rfs
KomaMRI.get_flip_angles
KomaMRI.get_ADC_on
KomaMRI.get_bvalue
KomaMRI.get_Bmatrix
KomaMRI.get_qvector
KomaMRI.get_M0_M1_M2
KomaMRI.get_max_grad
KomaMRI.get_RF_types
KomaMRI.get_kspace
KomaMRI.get_Mmatrix
KomaMRI.get_SRmatrix
KomaMRI.δ2N
KomaMRI.write_diff_fwf
KomaMRI.read_diff_fwf
```

### `Grad`
```@docs
rotx
roty
rotz
Grad
Grad(::Function, ::Real, ::Int64)
show(::IO, ::Grad)
getproperty(::Vector{Grad}, ::Symbol)
dur(::Grad)
```

### `RF`
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

RF
show(::IO, ::RF)
getproperty(::Vector{RF}, ::Symbol)
dur(::RF)
KomaMRI.RF_fun
KomaMRI.get_flip_angle
KomaMRI.get_RF_center
```

### `ADC`
```@docs
ADC
getproperty(::Vector{ADC}, ::Symbol)
KomaMRI.get_sample_times
KomaMRI.get_sample_phase_compensation
```

### `Delay`
```@docs
Delay
show(::IO, ::Delay)
+(::Sequence, ::Delay)
```

## [Pulseq.jl](@id pulseq)

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

### `read_seq`
```@docs
KomaMRI.read_seq
```

## [JEMRIS.jl](@id jemris)

### `read_phantom_jemris`
```@docs
read_phantom_jemris
```

## [ISMRMRD.jl](@id ismrmrd)

### `rawSignalToISMRMRD`
```@docs
rawSignalToISMRMRD
```

## [PulseDesigner.jl](@id pulse-designer)

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

### `simulate`
```@docs
simulate
```

### `simulate_slice_profile`
```@docs
simulate_slice_profile
```

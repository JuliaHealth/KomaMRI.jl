# Home

```@contents
Pages = ["index.md"]
Depth = 3
```

## [Dataflow Graph](@id dataflow-graph)

![](assets/dataflow.svg)

## [datatypes](@id datatypes)

### `Mag`
```@docs
Mag
```

### `Phantom`
```@docs
Phantom
```

### `Scanner`
```@docs
Scanner
```

### `Sequence`
```@docs
Sequence
KomaMRI.is_ADC_on
KomaMRI.is_GR_on
KomaMRI.is_RF_on
```

### `Grad`
```@docs
Grad
```

### `RF`
```@docs
Spinor
Rz
Ry
Rx
KomaMRI.Rg
KomaMRI.Rφ
Q

RF
KomaMRI.RF_fun
KomaMRI.get_flip_angle
KomaMRI.get_RF_center
```

### `ADC`
```@docs
ADC
KomaMRI.get_sample_times
KomaMRI.get_sample_phase_compensation
```

### `Delay`
```@docs
Delay
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

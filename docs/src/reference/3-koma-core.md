# KomaMRICore

```@meta
CurrentModule = KomaMRICore
```

## Simulation functions

```@docs
simulate
simulate_slice_profile
default_sim_params
```

## Simulation methods

```@docs
Bloch
BlochMagnusConst1
BlochMagnusLin2
BlochMagnusMid2
BlochMagnusLinComm2
BlochMagnusQuad2
BlochMagnusQuad4
BlochMagnusGL2
BlochMagnusGL4
BlochMagnusBGL4
BlochMagnusBGL6
```

## Receive weighting

```@docs
receive_weighted_dictionary
```

## GPU helper functions

```@docs
get_backend
print_devices
gpu
cpu
f32
f64
fbig
```

## Signal to `RawAquisitionData` (MRD)

```@docs
signal_to_raw_data
```

## `SpinRepresentationState`'s

```@docs
Mag
```

## `Spinor` rotation matrix (RF excitation)

```@docs
Spinor
Q
Un
Rx
Ry
Rz
```

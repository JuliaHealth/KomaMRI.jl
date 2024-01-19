# API Documentation

This page provides documentation for the modules, structs, functions, methods, and additional components available when importing the **KomaMRI** package. It serves as a valuable reference when using the **Julia REPL** directly and when creating custom **Julia** scripts. Be sure not to overlook the section [How to read the API docs](api.md#How-to-read-the-API-docs), which contains important information for understanding the general structure of docstrings. The following is the table of contents for the API Documentation:

```@raw html
<p align="center"><img width="100%" src="../assets/koma-schema-subdirs.svg"/></p>
```

```@contents
Pages = ["api.md"]
Depth = 3
```

## How to read the API docs

The API documentation includes predefined "template patterns" to assist users in understanding how to use modules, structs, functions, methods, and all the necessary aspects to make the most of what **KomaMRI** has to offer.

These documentation "template patterns" are based on the J[Julia Blue Style](https://github.com/invenia/BlueStyle)  documentation and other **GitHub** repositories that deal with MRI topics. However, some custom considerations were added to enhance understanding and provide a broader perspective.

When you encounter a docstring documentation, it will have the following structure:

!!! docstring "KomaMRI.component_name — Component"
    ```
    out1, out2, ... = component_name(arg1, arg2, ...; kw1, kw2, ...)
    ```

    This is a brief description of what **component_name** does.

    !!! note
        Here can be placed a note if it is regarded necessary.
    
    **Arguments**
    * `arg1`: (`::type`, `=value`, `[unit]`, opts: [`opt1`, `opt2`, ...]) the description for the arg1
    * ...
    
    **Keywords**
    * `kw1`: (`::type`, `=value`, `[unit]`, opts: [`opt1`, `opt2`, ...]) the description for the kw1
    * ...
    
    **Returns**
    * `out1`: (`::type`, `=value`, `[unit]`, opts: [`opt1`, `opt2`, ...]) the description for the out1
    * ...
    
    **References**
    * Sometimes it is a good idea to put some references or links
    * ...
    
    **Examples**
    ```julia-repl
    julia> arg1, arg2, valkw1, valkw2 = 3.5, "hello", 1, true

    julia> out1, out2 = component_name(arg1, arg2; kw1=valkw1, kw2=valkw2)
    ```

The preceding docstring block will always start with the way the component is called (outputs = component_name(inputs), followed by a brief description of what the component does. If necessary, a **note** block will be displayed. In general, the following subsections are optional: **Arguments**, **Keywords**, **Returns**, **References**, and **Examples**, but they will be provided as needed. These subsections are self-explanatory, making it intuitive to understand their purpose.

Please note that every subitem in the sections **Arguments**, **Keywords**, and **Returns** represents variables. They include practical information along with a description. The information enclosed in parentheses is optional but highly useful when provided.
* `::type`: is the suggested type for the variable. If the input variable is of type `::type`, there won't be any issues, but it's always possible to test other subtypes. If the variable is an output, it will be forced to the type `::type` whenever possible.
* `=value`: sometimes, for the inputs, a default value is defined if it is not assigned by the user.
* `[unit]`: this is the suggested physical unit of measure for the variable. Everything will be fine if you stick with these units of measure.
* opts: [`opt1`, `opt2`, ...]: sometimes, the input value can only be interpreted if it is one of the predefined values.

## Structs

```@meta
CurrentModule = KomaMRIBase
```

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


## Sequence-related Functions

### `dur`
```@docs
dur
```

### `rotx`
```@docs
rotx
```

### `roty`
```@docs
roty
```

### `rotz`
```@docs
rotz
```

### `get_kspace`
```@docs
get_kspace
```

### `get_M1`
```@docs
get_M1
```

### `get_M2`
```@docs
get_M2
```

### `get_flip_angles`
```@docs
get_flip_angles
```

### `is_RF_on`
```@docs
is_RF_on
```

### `is_GR_on`
```@docs
is_GR_on
```

### `is_ADC_on`
```@docs
is_ADC_on
```


## Prebuilt Sequence Building Blocks

```@meta
CurrentModule = KomaMRIBase
```

### `PulseDesigner`
```@docs
PulseDesigner
```

### `PulseDesigner.RF_hard`
```@docs
PulseDesigner.RF_hard
```

### `PulseDesigner.RF_sinc`
```@docs
PulseDesigner.RF_sinc
```

### `PulseDesigner.EPI`
```@docs
PulseDesigner.EPI
```

### `PulseDesigner.radial_base`
```@docs
PulseDesigner.radial_base
```

### `PulseDesigner.spiral_base`
```@docs
PulseDesigner.spiral_base
```

### `PulseDesigner.EPI_example`
```@docs
PulseDesigner.EPI_example
```


## Input/Output

```@meta
CurrentModule = KomaMRIFiles
```

### `read_seq`
```@docs
read_seq
```

### `read_phantom_jemris`
```@docs
read_phantom_jemris
```

### `read_phantom_MRiLab`
```@docs
read_phantom_MRiLab
```


## Functor Functions

```@meta
CurrentModule = KomaMRICore
```

**KomaMRICore** uses the Julia packages **Functors.jl** and **Adapt.jl**, which allows us to easily convert large nested structures. This could be used to send a custom type to the GPU (`gpu`), back to the CPU (`cpu`), change its floating point precision (`f32` and `f64`), etc. Four custom types in Koma implement these features: `Phantom`, `DiscreteSequence`, `Spinor`, and `Mag`. More information about these conversion operations is presented below.

### `f32`
```@docs
f32
```

### `f64`
```@docs
f64
```

### `cpu`
```@docs
cpu
```

### `gpu`
```@docs
gpu
```


## Simulation

### `default_sim_params`
```@docs
default_sim_params
```

### `simulate`
```@docs
simulate
```

### `simulate_slice_profile`
```@docs
simulate_slice_profile
```

### `signal_to_raw_data`
```@docs
signal_to_raw_data
```

### `Mag`
```@docs
Mag
```


## Plots

```@meta
CurrentModule = KomaMRIPlots
```

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

### `plot_M0`
```@docs
plot_M0
```

### `plot_M1`
```@docs
plot_M1
```

### `plot_M2`
```@docs
plot_M2
```

### `plot_eddy_currents`
```@docs
plot_eddy_currents
```

### `plot_slew_rate`
```@docs
plot_slew_rate
```

### `plot_signal`
```@docs
plot_signal
```

### `plot_image`
```@docs
plot_image
```

### `plot_dict`
```@docs
plot_dict
```

### `plot_seqd`
```@docs
plot_seqd
```


## UI

```@meta
CurrentModule = KomaMRI
```

### `KomaUI`
```@docs
KomaUI
```

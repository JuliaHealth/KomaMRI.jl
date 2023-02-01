
This page shows the documentation for the modules, structs, functions, methods and additional components available when importing the **KomaMRI.jl** package. It is very useful for reference when using directly the Julia REPL and when creating custom Julia scripts. Please, don't miss out the section [How to read the API docs](api.md#How-to-read-the-API-docs) which has important considerations to understand general aspects of the docstring structure. The following are the contents of the API Documentation:

```@contents
Pages = ["api.md"]
Depth = 3
```

## How to read the API docs

The API documentation has predefined "template patterns" which are meant to help the user to understand how to use modules, structs, functions, methods and every aspect necessary in order to take advantage of all the possibilities that **KomaMRI.jl** offers.

These documentation "template patterns" are based from the [Julia Blue Style](https://github.com/invenia/BlueStyle) documentation and other github repositories that works with MRI topics. However, some custom considerations were added for a better understanding and a wider perspective.

Whenever you see a docstring documentation, it will have the following structure:

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

The previous docstring block will always have first the way how it is called the component (outputs = component_name(inputs) and next a brief description about what the component does. Then a **note** block will be displayed if necessary. As general rule, the next subsections are optional: **Arguments**, **Keywords**, **Returns**, **References** and **Examples**, however they will be displayed whenever necessary. These subsections are self-explanatory, so it is intuitive to figure out what are they meant for.

Note that every subitem in the sections **Arguments**, **Keywords** and **Returns** are variables. They have practical information enclosed in parentheses plus a description. They information in parentheses is optional but very useful if it is present:
* `::type`: is the suggested type of the variable. If the input variable is of type `::type`, then nothing can go wrong, but it is always possible to test other types. If the variable is an output, then it will always try to be forced to the type `::type`.
* `=value`: sometimes for the inputs the is defined a default value if it is not assigned by the user.
* `[unit]`: this is the suggested physical unit of measure of the variable. Everything it is going to be fine if you are stick with these units of measure.
* opts: [`opt1`, `opt2`, ...]: sometimes the input value can only be interpreted if it is part of some predefined values.

## Structs

```@meta
CurrentModule = KomaMRICore
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

## Read Data

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

### `signal_to_raw_data`
```@docs
signal_to_raw_data
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

### `plot_M0`
```@docs
plot_M0
```

### `plot_signal`
```@docs
plot_signal
```

### `plot_image`
```@docs
plot_image
```

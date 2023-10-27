# Sequence Definition

This section dives into some details about how a sequence is constructed. The sequence definition of **KomaMRI** is strongly realated with the [Pulseq](https://pulseq.github.io/index.html) definition. After you read this section, you should be able to construct your own **Sequence** structs to perform your custom simulations using the **KomaMRI** package.

## Sequence Overview

Let's introduce the following simple sequence figure to extend the ideas from a visual example to a more general sequence definition:

```@raw html
<p align="center"><img width="80%" src="../assets/sequence-diagram.svg"/></p>
```

A **sequence** can be thought as and ordered concatenation of blocks over time. A block is a **sequence** of length 1. Every block is composed by an **RF** pulse, the ``(x,y,z)`` **gradients**,  and the **acquisition** of the samples. There is also a time **duration** associated to each block. For short, we are going to refer to these components like so:

```math
\begin{matrix*}[l]
seq[i]      &: & \text{block $i$ of the sequence} \\
seq.RF[i]   &: & \text{RF pulse at the $i$ block} \\
seq.GR.x[i] &: & \text{gradient x at the $i$ block} \\
seq.GR.y[i] &: & \text{gradient y at the $i$ block} \\
seq.GR.z[i] &: & \text{gradient z at the $i$ block} \\
seq.ADC[i]  &: & \text{acquisition at the $i$ block} \\
seq.DUR[i]  &: & \text{duration at the $i$ block}
\end{matrix*}
```

The best way to understand the **Sequence** struct of **KomaMRI** is by directly seeing the source code where this struct is defined:
```julia
mutable struct Sequence
    GR::Array{Grad,2}      # Sequence in (X, Y and Z) and time
    RF::Array{RF,2}        # RF pulses in coil and time
    ADC::Array{ADC,1}      # ADC in time
    DUR::Vector            # Duration of each block
    DEF::Dict{String,Any}  # Dictionary with information relevant to the reconstructor
end
```

As you can see, a **Sequence** struct contains 5 field names: ''DEF'' contains information for reconstruction steps (so it is not mandatory to fill it), ''DUR'' is a vector that contains the time durations of each block, ''ADC'' is also a vector with the acquisition samples for every block (an vector of **ADC** structs), ''GR'' is a 2D matrix which 3 rows representing the x-y-z gradients and columns having the samples of each block (a matrix of **Grad** structs) and ''RF'' is also a 2D matrix where each row represents a different coil and the columns are for different block samples too (a matrix of **RF** structs). The **RF**, **Grad** and **ADC** are MRI events that will be explained in the section [Events Definitions](events.md).

!!! warning
    So far, **KomaMRI** only can manage 1 coil for RF excitations

In order to understand how a **Sequence** struct can be manipulated in **Julia**, let's use the EPI sequence example. You can display basic information of the **Sequence** variable in the **Julia REPL**:
```julia-repl
julia> seq = PulseDesigner.EPI_example()
Sequence[ τ = 62.846 ms | blocks: 204 | ADC: 101 | GR: 205 | RF: 1 | DEF: 5 ]
```

As you can see, this **Sequence** has 204 blocks, 1 of these blocks has an **RF** struct with values different from zero, there are 205 number of **Grad** structs considering the x-y-z components, 101 **ADC** structs acquire samples of some blocks and 62.846 ms is the total time duration of the complete **Sequence**.

To display the sequence in an graph, we can use the **plot\_seq()** function:
```julia-repl
julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="100%" src="../assets/seq-epi-example-full.svg"/></p>
```

In this way you can see exactly where are present the **RF**, **Grad** and **ADC** structs.

You can filter the information for the **RF**, **Grad**, **ADC** and **DUR** field names of a **Sequence** by simply accessing to them with the ''dot'' notation, so you can display helpful information of how the **Sequence** struct is organized:
```julia-repl
julia> seq.RF
1×204 Matrix{RF}:
 ⊓(0.5872 ms)  ⇿(0.0 ms)  ⇿(0.0 ms)  …  ⇿(0.0 ms)  ⇿(0.0 ms)   

julia> seq.GR
3×204 Matrix{Grad}:
 ⇿(0.5872 ms)  ⊓(0.4042 ms)  ⊓(0.4042 ms)  …  ⇿(0.2062 ms)  ⊓(0.4042 ms)  ⊓(0.4042 ms)
 ⇿(0.5872 ms)  ⊓(0.4042 ms)  ⇿(0.4042 ms)     ⊓(0.2062 ms)  ⇿(0.4042 ms)  ⊓(0.4042 ms)
 ⇿(0.5872 ms)  ⇿(0.0 ms)     ⇿(0.0 ms)        ⇿(0.0 ms)     ⇿(0.0 ms)     ⇿(0.0 ms)

julia> seq.ADC
204-element Vector{ADC}:
 ADC(0, 0.0, 0.0, 0.0, 0.0)
 ADC(0, 0.0, 0.0, 0.0, 0.0)
 ADC(101, 0.00019999999999999998, 0.00010211565434713023, 0.0, 0.0)
 ⋮
 ADC(101, 0.00019999999999999998, 0.00010211565434713023, 0.0, 0.0)
 ADC(0, 0.0, 0.0, 0.0, 0.0)

julia> seq.DUR
204-element Vector{Float64}:
 0.0005871650124959989
 0.0004042313086942605
 0.0004042313086942605
 ⋮
 0.0004042313086942605
 0.0004042313086942605
```

Additionally, you can access to a subset of blocks in a **Sequence** by slicing or indexing, in which case the result will also be a **Sequence** struct and thus you can perform the same normal operations for a **Sequence**. For example, if you are interested in analyzing just the first 11 blocks, you can make something like this:
```julia-repl
julia> seq[1:11]
Sequence[ τ = 3.837 ms | blocks: 11 | ADC: 5 | GR: 11 | RF: 1 | DEF: 5 ]

julia> seq[1:11].GR
3×11 Matrix{Grad}:
 ⇿(0.5872 ms)  ⊓(0.4042 ms)  ⊓(0.4042 ms)   …  ⊓(0.4042 ms)  ⇿(0.2062 ms)  ⊓(0.4042 ms)
 ⇿(0.5872 ms)  ⊓(0.4042 ms)  ⇿(0.4042 ms)      ⇿(0.4042 ms)  ⊓(0.2062 ms)  ⇿(0.4042 ms)
 ⇿(0.5872 ms)  ⇿(0.0 ms)     ⇿(0.0 ms)        ⇿(0.0 ms)     ⇿(0.0 ms)     ⇿(0.0 ms)

julia> plot_seq(seq[1:11])
```
```@raw html
<p align="center"><img width="100%" src="../assets/seq-epi-example-some-blocks.svg"/></p>
```

## Concatenation of Sequences 

We can concatenate sequences together side by side. The example bellow shows how to concatenate sequences:
```julia-repl
julia> s = PulseDesigner.EPI_example()[1:11]
Sequence[ τ = 3.837 ms | blocks: 11 | ADC: 5 | GR: 11 | RF: 1 | DEF: 5 ]

julia> seq = s + s + s
Sequence[ τ = 11.512 ms | blocks: 33 | ADC: 15 | GR: 33 | RF: 3 | DEF: 5 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/seq-concatenation.svg"/></p>
```
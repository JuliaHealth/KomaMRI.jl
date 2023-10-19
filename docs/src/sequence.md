# Sequence Definition

This section dives into some details about how a sequence is constructed. After you read this section, you should be able to construct your own **Sequence** structs to perform your custom simulations using the **KomaMRI** package.

## Sequence

Let's introduce the following simple sequence figure to extend the ideas from a visual example to a more general sequence definition:

```@raw html
<p align="center"><img width="60%" src="../assets/sequence-diagram.svg"/></p>
```

A **sequence** can be thought as and ordered concatenation of blocks over time. Every block is composed by an **RF** pulse, the ``(x,y,z)`` **gradients**,  and the **acquisition** of the samples. There is also a time **duration** associated to each block. For short, we are going to refer to these components like so:

```math
\begin{matrix*}[l]
i          &: & \text{sequence block ID} \\
RF[i]      &: & \text{RF pulse at the $i$ block} \\
G_j[i]     &: & \text{gradients at the $i$ block}, \: \forall j \in \{x,y,z\} \\
ADC[i]     &: & \text{acquisition at the $i$ block} \\
DUR[i]     &: & \text{duration at the $i$ block}
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

As you can see, a **Sequence** struct contains 5 field names: ''DEF'' contains information for reconstruction steps (so it is not mandatory to fill it), ''DUR'' is a vector that contains the time durations of each block, ''ADC'' is also a vector with the acquisition samples for every block (an vector of **ADC** structs), ''GR'' is a 2D matrix which 3 rows representing the x-y-z gradients and columns having the samples of each block (a matrix of **Grad** structs) and ''RF'' is also a 2D matrix where each row represents a different coil and the columns are for different block samples too (a matrix of **RF** structs). The **RF**, **Grad** and **ADC** structs will be explained later in their corresponding subsections.

!!! note
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

Additionally, you can access to limited amount of blocks in a **Sequence** by considering block indices or ranges, in which case the result will also be a **Sequence** struct and thus you can perform the same normal operations of a **Sequence**. For example, if you are interested in analyzing just the first 11 blocks, you can make something like this:
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


## Event Waveforms

We refer to **RF**, **Grad** and **ADC** as ''events''. As we already know, a **Sequence** struct contains field names that access to arrays of **RF**, **Grad** and **ADC** structs. So, in order to create a **Sequence** you need to understand how to create these basic events. In the image bellow, we show a summary of how you can define **RF**, **Grad** and **ADC** events:

```@raw html
<p align="center"><img width="90%" src="../assets/event-shapes-rf-zoh.svg"/></p>
```
```@raw html
<!--
<p align="center"><img width="90%" src="../assets/event-shapes.svg"/></p>
-->
```

As a general overview, you can see that the **RF** and **Grad** structs has 3 ways of defining waveforms: a simple pulse/trapezoidal, evenly time spaced waveform (uniformly-sampled) and non-uniform time sampled waveform (time-shaped); whereas the **ADC** just have one way to define its struct. It's important to notice that the way the samples of the **RF** and **Grad** structs are interpolated by **KomaMRI** are different; on one hand for **RF** is used a zero-order-hold sampling, on the other hand for **Grads** is used linear interpolation. 

!!! note
    The way in which the **RF** interpolation is interpreted by **KomaMRI** will change in future versions considering linear interpolation between 2 consecutive samples instead of the present zero-order-hold, just like is being done with the **Grad** struct.

In the following subsections, we are going to explain more in detail the parameters to define events and how to create a **Sequence** with the **RF**, **Grad** and **ADC** events.


## RF

The **RF** struct is defined as follows in the source code of **KomaMRI**:
```julia
mutable struct RF
    A
    T
    Δf
    delay::Float64
    ...
end
```

As you can see, it has 4 field names: ''A'' defines amplitude, ''T'' defines duration time, ''delay'' is the distance between the 0 time and the first waveform sample and ''Δf'' is the displacement respect to the main field carrier frequency (this is for advanced users).

''A'' and ''T'' can be numbers or vectors of numbers. Depending on the length of the ''A'' and ''T'', **KomaMRI** interprets different waveforms: 
* Pulse Waveform: A and T are numbers
* Uniformly-Sampled Waveform: A is a vector and T is a number
* Time-Shaped Waveform: A and T are both vectors with the same length (zero-order-hold)

Let's see some basic examples on how to create these types of **RF** structs and put them inside a **Sequence** struct. The examples should be self explanatory.

### RF Pulse Waveform

```julia-repl
julia> A, T, delay =  10e-3, 0.5e-3, 0.1e-3;

julia> rf = RF(A, T, 0, delay)
←0.1 ms→ RF(10000.0 uT, 0.5 ms, 0.0 Hz)

julia> seq = Sequence(); seq += rf; seq = seq[2:end]
Sequence[ τ = 0.6 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-rf-pulse-waveform.svg"/></p>
```

### RF Uniformly-Sampled Waveform

```julia-repl
julia> A, T, delay =  10e-3*[0; -0.1; 0.2; -0.5; 1; -0.5; 0.2; -0.1; 0], 0.5e-3, 0.1e-3;

julia> rf = RF(A, T, 0, delay)
←0.1 ms→ RF(∿ uT, 0.5 ms, 0.0 Hz)

julia> seq = Sequence(); seq += rf; seq = seq[2:end]
Sequence[ τ = 0.6 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-rf-uniformly-samples-waveform.svg"/></p>
```

### RF Time-Shaped Waveform

```julia-repl
julia> A =  10e-3*[0; -0.1; 0.2; -0.5; 1; -0.5; 0.2; -0.1; 0];

julia> T, delay = 0.5e-3*[0.1; 0.2; 0.1; 0.2; 0.2; 0.1; 0.2; 0.1; 0.15], 0.1e-3;

julia> rf = RF(A, T, 0, delay)
←0.1 ms→ RF(∿ uT, 0.675 ms, 0.0 Hz)

julia> seq = Sequence(); seq += rf; seq = seq[2:end]
Sequence[ τ = 0.775 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-rf-time-shaped-waveform.svg"/></p>
```


## Gradient

The **Grad** struct is defined as follows in the source code of **KomaMRI**:
```julia
mutable struct Grad
    A
    T
    rise::Real
    fall::Real
    delay::Real
    ...
end
```

As you can see, it has 5 field names: ''A'' defines amplitude, ''T'' defines duration time, ''delay'' is the distance between the 0 time and the first waveform sample, ''rise'' and ''fall'' are the time durations of the first and last gradient ramps.

Just like the **RF**, ''A'' and ''T'' in the **Grad** struct can be numbers or vectors of numbers. Depending on the length of the ''A'' and ''T'', **KomaMRI** interprets different waveforms: 
* Trapezoidal Waveform: A and T are numbers
* Uniformly-Sampled Waveform: A is a vector and T is a number
* Time-Shaped Waveform: A and T are both vectors, A has one sample more the T (linear interpolation)

Let's see some basic examples on how to create these types of **Grad** structs and put them inside a **Sequence** struct considering just the ''x'' component of the gradients. The examples should be self explanatory.

### Gradient Trapezoidal Waveform

```julia-repl
julia> A, T, delay, rise, fall =  50*10e-6, 5e-3, 2e-3, 1e-3, 1e-3;

julia> gr = Grad(A, T, rise, fall, delay)
←2.0 ms→ Grad(0.5 mT, 0.5 ms, ↑1.0 ms, ↓1.0 ms)

julia> seq = Sequence([gr])
Sequence[ τ = 9.0 ms | blocks: 1 | ADC: 0 | GR: 1 | RF: 0 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-gr-trapezoidal-waveform.svg"/></p>
```

### Gradient Uniformly-Sampled Waveform

```julia-repl
julia> A = 50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1];

julia> T, delay, rise, fall = 5e-3, 2e-3, 0, 1e-3;

julia> gr = Grad(A, T, rise, fall, delay)
←2.0 ms→ Grad(∿ mT, 5.0 ms, ↑0.0 ms, ↓1.0 ms)

julia> seq = Sequence([gr])
Sequence[ τ = 8.0 ms | blocks: 1 | ADC: 0 | GR: 1 | RF: 0 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-gr-uniformly-sampled-waveform.svg"/></p>
```

### Gradient Time-Shaped Waveform

```julia-repl
julia> A = 50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1];

julia> T = 5e-3*[0.2; 0.1; 0.3; 0.2; 0.1; 0.2; 0.3; 0.2; 0.1];

julia> delay, rise, fall = 2e-3, 1e-3, 1e-3;

julia> gr = Grad(A, T, rise, fall, delay)
←2.0 ms→ Grad(∿ mT, 7.75 ms, ↑0.0 ms, ↓1.0 ms)

julia> seq = Sequence([gr])
Sequence[ τ = 10.75 ms | blocks: 1 | ADC: 0 | GR: 1 | RF: 0 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-gr-time-shaped-waveform.svg"/></p>
```

## ADC

The **ADC** struct is defined as follows in the source code of **KomaMRI**:
```julia
mutable struct ADC
    N::Int64
    T::Float64
    delay::Float64
    Δf::Float64
    ϕ::Float64
    ...
end
```

As you can see, it has 5 field names: ''N'' defines number of samples, ''T'' defines total acquisition duration, ''delay'' is the distance between the 0 time and the first sampled signal, ''Δf'' and ''ϕ' are factor to correct signal acquisition (for advanced users).

Let's see a basic example on how to define an **ADC** struct and put it inside a **Sequence** struct:
```julia-repl
julia> N, T, delay =  16, 5e-3, 1e-3;

julia> adc = ADC(N, T, delay)
ADC(16, 0.005, 0.001, 0.0, 0.0)

julia> seq = Sequence(); seq += adc; seq = seq[2:end]
Sequence[ τ = 6.0 ms | blocks: 1 | ADC: 1 | GR: 0 | RF: 0 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-adc.svg"/></p>
```

## Event Combination

We call ''event combination'' to the action of putting inside a sequence of one block multiple events. The example bellow shows how to combine an **RF** struct, 3 **Grad** structs for the x-y-z components and an **ADC** struct in a sequence of one block: 
```julia
# Define an RF struct
A, T, delay =  1e-6*[0; -0.1; 0.2; -0.5; 1; -0.5; 0.2; -0.1; 0], 0.5e-3, 0.1e-3;
rf = RF(A, T, 0, delay)

# Define a Grad struct for Gx
A, T, delay, rise, fall =  50*10e-6, 5e-3, 2e-3, 1e-3, 1e-3
gx = Grad(A, T, rise, fall, delay)

# Define a Grad struct for Gy
A = 50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1]
T, delay, rise, fall = 5e-3, 2e-3, 0, 1e-3;
gy = Grad(A, T, rise, fall, delay)

# Define a Grad struct for Gz
A = 50*10e-6*[0; 0.5; 0.9; 1; 0.9; 0.5; 0; -0.5; -0.9; -1]
T = 5e-3*[0.2; 0.1; 0.3; 0.2; 0.1; 0.2; 0.3; 0.2; 0.1]
delay, rise, fall = 2e-3, 1e-3, 1e-3
gz = Grad(A, T, rise, fall, delay)

# Define an ADC struct
N, T, delay =  16, 5e-3, 1e-3
adc = ADC(N, T, delay)
```
```julia-repl
julia> seq = Sequence([gx; gy; gz;;], [rf;;], [adc])
Sequence[ τ = 12.5 ms | blocks: 1 | ADC: 1 | GR: 3 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-combination.svg"/></p>
```

It is important to note that once the struct events are defined, in order to create a sequence of one block, is necessary to put 2D matrices of **Grad** and **RF** structs and a vector of **ADC** structs as arguments in the **Sequence()** constructor.

## Sequence Concatenation

We call ''sequence concatenation'' to the action of putting sequences together side by side. The example bellow shows how to concatenate sequences:
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
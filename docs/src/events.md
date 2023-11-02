# Event Definitions

We refer to **RF**, **Grad** and **ADC** as ''events''. This section deals with some details about how events can be defined and how they can be manipulated inside a **Sequence** struct.

## Events Overview

As we already know, a **Sequence** struct contains field names that access to arrays of **RF**, **Grad** and **ADC** structs. So, in order to create a **Sequence** you need to understand how to create these basic events. In the image bellow, we show a summary of how you can define **RF**, **Grad** and **ADC** events:

```@raw html
<p align="center"><img width="90%" src="../assets/event-shapes.svg"/></p>
```

As a general overview, you can see that the **RF** and **Grad** structs has 3 ways of defining waveforms: a simple pulse/trapezoidal, evenly time spaced waveform (uniformly-sampled) and non-uniform time sampled waveform (time-shaped); whereas the **ADC** just have one way to define its struct. It's important to notice that the way the samples of the **RF** and **Grad** structs are interpolated by **KomaMRI** are different; on one hand for **RF** is used a zero-order-hold sampling, on the other hand for **Grads** is used linear interpolation. 

!!! warning
    The way in which the **RF** interpolation is interpreted by **KomaMRI** will change in future versions considering linear interpolation between 2 consecutive samples instead of the present zero-order-hold, just like is being done with the **Grad** struct.

In the following subsections, we are going to explain more in detail the parameters to define events and how to create a **Sequence** with the **RF**, **Grad** and **ADC** events.


## RF

The **RF** struct is defined as follows in the source code of **KomaMRI**:
```julia
mutable struct RF
    A
    T
    Δf
    delay::Real
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
julia> tl = -3:0.2:-0.2; tr = 0.2:0.2:3;

julia> A = (10e-3)*[sin.(π*tl)./(π*tl); 1; sin.(π*tr)./(π*tr)];

julia> T, delay = 0.5e-3, 0.1e-3;

julia> rf = RF(A, T, 0, delay)
←0.1 ms→ RF(∿ uT, 0.5 ms, 0.0 Hz)

julia> seq = Sequence(); seq += rf; seq = seq[2:end]
Sequence[ τ = 0.6 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-rf-uniformly-sampled-waveform.svg"/></p>
```

### RF Time-Shaped Waveform

```julia-repl
julia> tl = -4:0.2:-0.2; tr = 0.2:0.2:4

julia> A = (10e-3)*[sin.(π*tl)./(π*tl); 1; sin.(π*tr)./(π*tr)]

julia> T = [0.05e-3*ones(length(tl)); 2e-3; 0.05e-3*ones(length(tl))]

julia> delay = 0.1e-3;

julia> T, delay = 0.5e-3*[0.1; 0.2; 0.1; 0.2; 0.2; 0.1; 0.2; 0.1; 0.15], 0.1e-3;

julia> rf = RF(A, T, 0, delay)
←0.1 ms→ RF(∿ uT, 4.0 ms, 0.0 Hz)

julia> seq = Sequence(); seq += rf; seq = seq[2:end]
Sequence[ τ = 4.1 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

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
julia> t = 0:0.25:7.5

julia> A = 10*10e-6 * sqrt.(π*t) .* sin.(π*t)

julia> T = 10e-3;

julia> delay, rise, fall = 1e-3, 0, 1e-3;

julia> gr = Grad(A, T, rise, fall, delay)
←1.0 ms→ Grad(∿ mT, 10.0 ms, ↑0.0 ms, ↓1.0 ms)

julia> seq = Sequence([gr])
Sequence[ τ = 12.0 ms | blocks: 1 | ADC: 0 | GR: 1 | RF: 0 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-gr-uniformly-sampled-waveform.svg"/></p>
```

### Gradient Time-Shaped Waveform

```julia-repl
julia> A = 50*10e-6*[1; 1; 0.8; 0.8; 1; 1];

julia> T = 1e-3*[5; 0.2; 5; 0.2; 5];

julia> delay, rise, fall = 1e-3, 1e-3, 1e-3;

julia> gr = Grad(A, T, rise, fall, delay)
←1.0 ms→ Grad(∿ mT, 15.4 ms, ↑1.0 ms, ↓1.0 ms)

julia> seq = Sequence([gr])
Sequence[ τ = 10.75 ms | blocks: 1 | ADC: 0 | GR: 1 | RF: 0 | DEF: 0 ]

julia> plot_seq(seq)
```
```@raw html
<p align="center"><img width="90%" src="../assets/event-gr-time-shaped-waveform.svg"/></p>
```

## ADC

The **ADC** struct is defined in the **KomaMRI** source code as follows:
```julia
mutable struct ADC
    N::Integer
    T::Real
    delay::Real
    Δf::Real
    ϕ::Real
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

!!! warning
    The **Pulseq** definition of the ADC is different from the **KomaMRI** definition. **KomaMRI** considers ADC samples at `delay` time, at `delay + T`, and with equispaced positions for remaining samples (except for and ADC with 1 sample in which is at the position `delay + T/2`). Meanwhile **Pulseq** regards equispaced samples without considering the `delay` nor the `delay + T` points.    

## Combination of Events

We can insert multiple events inside a sequence of one block. The example bellow shows how to combine an **RF** struct, 3 **Grad** structs for the x-y-z components and an **ADC** struct in a sequence of one block: 
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

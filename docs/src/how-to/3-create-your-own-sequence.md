# Create Your Own Sequence

!!! warning
    This section is currently under construction, and some details on how to construct a Sequence may be missing.

This is an example of how to create a **Sequence** struct:

```julia
# Export necessary modules
using KomaMRI

# Create the function that creates a phantom
function sequence_example(FOV::Real, N::Integer)

    # Define initial paramters (TODO: consider when N is even)
    sys = Scanner()
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	Ga ≥ sys.Gmax ? error("Ga=$(Ga*1e3) mT/m exceeds Gmax=$(Gmax*1e3) mT/m, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	ϵ1 = Δτ/(Δτ+ζ)

	# EPI base
	EPI = Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2],  #Gx
	 	[mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2])) #Gy
	EPI.ADC = [mod(i,2)==1 ? ADC(0,Δτ,ζ) : ADC(N,Ta,ζ) for i=0:2*Ny-2]

	# Pre-wind and wind gradients
	ϵ2 = Ta/(Ta+ζ)
    PHASE =   Sequence(reshape(1/2*[Grad(      -Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) # This needs to be calculated differently
	DEPHASE = Sequence(reshape(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) # for even N
	seq = PHASE + EPI + DEPHASE

	# Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")

    # Return the sequence
	return seq
end

# Call the function to create a sequence
FOV, N = 23e-2, 101
seq = sequence_example(FOV, N)

# Plot the sequence in time and its kspace
plot_seq(seq; range=[0 30])
plot_kspace(seq)
```
```@raw html
<object type="text/html" data="../assets/create-your-own-sequence-time.html" style="width:50%; height:420px;"></object><object type="text/html" data="../assets/create-your-own-sequence-kspace.html" style="width:50%; height:420px;"></object>
```

## Block-Oriented Construction

Koma sequences can also be built in the same style as MATLAB Pulseq or PyPulseq:
define reusable events, then add one physical block at a time. This is often the
clearest way to write pulse programs with loops.

:::tabs

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz)
seq.addBlock(gx, adc)
```

== PyPulseq

```python
seq.add_block(rf, gz)
seq.add_block(gx, adc)
```

== KomaMRI

```julia
addblock!(seq, rf; z=gz)
addblock!(seq, adc; x=gx)
```

:::

Pulseq gradient events carry their channel. Koma `Grad` events are axis-neutral,
so the axis is chosen when the block is added.

```julia
seq = Sequence()

for ky in 1:Ny
    addblock!(seq, rf; z=gz)
    addblock!(seq, adc; x=gx, y=phase_blip(ky))
end
```

The axis keywords `x=`, `y=`, and `z=` place gradients in the corresponding
gradient row. RF, ADC, and extensions such as `LabelSet`, `LabelInc`, and
`Trigger` are positional block events. `Delay(T)` can be used as a standalone
delay block, or inside a block to set the minimum block duration. `Duration(T)`
sets the exact block duration and errors if any event is longer.

```julia
addblock!(seq, rf, LabelSet(ky, "LIN"), Duration(TR); z=gz)
```

This mirrors Pulseq's two duration styles: `Delay(T)` is like
`mr.makeDelay(T)`, while `Duration(T)` is like MATLAB Pulseq's numeric
`seq.addBlock(T, ...)` argument.

For longer sequences, `@addblocks` keeps the Pulseq mental model while making
composition concise:

!!! warning
    Do not use plain `seq += ...` in long loops without `@addblocks`. Outside the
    macro, `+=` is just `seq = seq + ...`, so it copies the accumulated sequence
    on every iteration and becomes very inefficient for large pulse programs.

```julia
seq = Sequence()

@addblocks for ky in 1:Ny
    seq += (rf, LabelSet(ky, "LIN"), Duration(TR), z=gz)
    seq += (adc, x=gx, y=phase_blip(ky))
end
```

Inside `@addblocks`, a tuple means "one block" and `+` means "append these
chunks in order". Incoming events and sequence chunks are copied when appended,
so reusing `rf`, `gz`, or a readout chunk in a loop does not connect old blocks
to later edits.

Variable block contents can be spliced, matching the MATLAB `contents{:}` style:

```julia
contents = (Duration(TR), LabelInc(1, "ECO"))
@addblocks seq += (contents..., x=gx, y=gy, z=gz)
```

```julia
bSSFP = Sequence()

@addblocks begin
    for (i, f) in enumerate(ramp_factors)
        bSSFP += ((-1)^(i - 1) * f) * rf + ramp_readout
    end

    for ky in 1:Ny
        bSSFP += (-1)^(ky - 1) * rf + readout(ky)
    end
end
```

Use `@addblock` for one statement and `@addblocks` for a loop or `begin ... end`
block.

:::tabs

== Source

```julia
@addblock seq += (rf, LabelSet(line, "LIN"), Duration(TR), z=gz)
```

== Equivalent explicit code

```julia
addblock!(seq, rf, LabelSet(line, "LIN"), Duration(TR); z=gz)
```

:::

:::tabs

== Source

```julia
@addblocks for ky in 1:Ny
    seq += rf_prep + readout(ky)
    counter += 1
end
```

== Equivalent explicit code

```julia
for ky in 1:Ny
    append!(seq, rf_prep)
    append!(seq, readout(ky))
    counter = counter + 1
end
```

:::

The macro is only syntax for explicit block appends. `Sequence` left-hand sides
append blocks in place, while non-`Sequence` `+=` expressions keep normal Julia
meaning.

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
@addblock seq = (rf, z=gz)
@addblock seq += (adc, x=gx)
```

:::

Koma gradients are axis-neutral. Choose the axis when adding the block.

```julia
seq = Sequence()

for ky in 1:Ny
    @addblock seq += (rf, z=gz)
    @addblock seq += (adc, x=gx, y=phase_blip(ky))
end
```

RF, ADC, and extensions are positional. Gradients use `x=`, `y=`, or `z=`.

```julia
@addblock seq += (rf, LabelSet(ky, "LIN"), Duration(TR), z=gz)
```

Use `Delay(T)` for a minimum duration. Use `Duration(T)` for an exact duration.

```julia
@addblock seq += (rf, Delay(TR), z=gz)     # at least TR
@addblock seq += (rf, Duration(TR), z=gz)  # exactly TR, or errors
```

Use `@addblocks` in loops.

!!! warning
    Plain `seq += chunk` in long loops copies the accumulated sequence each
    iteration. In a 1,000-block loop, this can be more than 800x slower than
    `@addblocks`.

```julia
# Slow
for ky in 1:Ny
    seq += readout(ky)
end

# Fast
@addblocks for ky in 1:Ny
    seq += readout(ky)
end
```

```julia
seq = Sequence()

@addblocks for ky in 1:Ny
    seq += (rf, LabelSet(ky, "LIN"), Duration(TR), z=gz)
    seq += (adc, x=gx, y=phase_blip(ky))
end
```

Inside `@addblocks`, a tuple is one block. `+` appends chunks in order.

Splice variable block contents:

```julia
contents = (Duration(TR), LabelInc(1, "ECO"))
@addblocks seq += (contents..., x=gx, y=gy, z=gz)
```

Group gradient axes:

```julia
rf_event = RF(rf_waveform, rf_intervals, Δf, rf_delay)
G_ss = (x=Grad(Gx, T, ζ), y=Grad(Gy, T, ζ), z=Grad(Gz, T, ζ))
G_rew = (x=Grad(Gx_rew, T_rew, ζ), y=Grad(Gy_rew, T_rew, ζ), z=Grad(Gz_rew, T_rew, ζ))

@addblock excitation = (rf_event; G_ss...) + (; G_rew...)
```

Name reusable chunks:

```julia
@addblock RO = (ADC(N_ro, T_adc, ζ), x=Grad(G_ro, T_ro, ζ))
@addblock PRE = (x=Grad(Gx_pre, T_pre, ζ_pre), y=Grad(Gy_pre, T_pre, ζ_pre))
@addblock line = PRE + RO
```

Appending copies incoming events and chunks, so reused chunks do not share mutable
events with appended blocks.

```julia
bSSFP = Sequence()

@addblocks for (i, f) in enumerate(ramp_factors)
    phase = cis(π * (i - 1))
    bSSFP += rf_excitation(f * phase) + phase * ramp_readout
end

@addblocks for ky in 1:Ny
    phase = cis(π * (length(ramp_factors) + ky - 1))
    bSSFP += rf_excitation(phase) + phase * readout(ky)
end
```

Complex scaling phase-shifts RF and ADC and leaves gradients unchanged.

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

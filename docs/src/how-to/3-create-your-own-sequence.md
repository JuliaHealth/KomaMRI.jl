# Build Sequences with `@addblock`

A Koma `Sequence` is a list of sequence blocks. Each block stores RF events,
one gradient per axis, ADC events, a duration, and extensions such as labels,
triggers, or rotations.

This page shows the block-oriented style for building pulse programs. It maps
closely to MATLAB Pulseq and PyPulseq `addBlock` / `add_block` code, but Koma can
also append named `Sequence` parts, such as a readout or prephaser, in one
statement.

Use `@addblock` to append block expressions:

:::tabs

== KomaMRI

```julia
seq = Sequence()  # or Sequence(sys) for scanner export checks
@addblock seq += (rf, z=gz)
```

== MATLAB Pulseq

```matlab
seq = mr.Sequence();
seq.addBlock(rf, gz)
```

== PyPulseq

```python
seq = pp.Sequence()
seq.add_block(rf, gz)
```

:::

After `+=`, each tuple becomes a `Sequence` of length one and its events are
copied into `seq`. Each `Sequence` on the right-hand side contributes all of its
blocks. `+` appends them left-to-right:

:::tabs

== KomaMRI

```julia
@addblock seq += (rf, z=gz) + (x=gx, adc) + (Delay(TR), LabelInc(1, "ECO"))
```

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz)
seq.addBlock(gx, adc)
seq.addBlock(mr.makeDelay(TR), mr.makeLabel('INC', 'ECO', 1))
```

== PyPulseq

```python
seq.add_block(rf, gz)
seq.add_block(gx, adc)
seq.add_block(pp.make_delay(TR), pp.make_label('ECO', 'INC', 1))
```

:::

## Gradient Axes

Koma `Grad` events are axis-neutral. Choose the axis when adding the block with
`x=`, `y=`, or `z=`. RF, ADC, and extensions are positional.

:::tabs

== KomaMRI

```julia
@addblock seq += (rf, LabelSet(ky, "LIN"), z=gz)
```

== MATLAB Pulseq

```matlab
seq.addBlock(rf, mr.makeLabel('SET', 'LIN', ky), gz)
```

== PyPulseq

```python
seq.add_block(rf, pp.make_label('LIN', 'SET', ky), gz)
```

:::

Named tuples are useful when a block has gradients on several axes:

:::tabs

== KomaMRI

```julia
slice_select = (x=gx_slice, y=gy_slice, z=gz_slice)
rewinder = (x=gx_rewind, y=gy_rewind, z=gz_rewind)

excitation = Sequence()
@addblock excitation += (rf; slice_select...) + (; rewinder...)
```

== MATLAB Pulseq

```matlab
sliceSelect = {gxSlice, gySlice, gzSlice};
rewinder = {gxRewind, gyRewind, gzRewind};

excitation = mr.Sequence();
excitation.addBlock(rf, sliceSelect{:})
excitation.addBlock(rewinder{:})
```

== PyPulseq

```python
slice_select = [gx_slice, gy_slice, gz_slice]
rewinder = [gx_rewind, gy_rewind, gz_rewind]

excitation = pp.Sequence()
excitation.add_block(rf, *slice_select)
excitation.add_block(*rewinder)
```

:::

## Block Duration

Use `Delay(T)` for a minimum block duration. Use `Duration(T)` for an exact block
duration; it errors if any RF, gradient, or ADC event is longer than `T`.

:::tabs

== KomaMRI

```julia
@addblock seq += (rf, Delay(TR), z=gz)     # at least TR
@addblock seq += (rf, Duration(TR), z=gz)  # exactly TR, or errors
```

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz, mr.makeDelay(TR))
seq.addBlock(TR, rf, gz)
```

== PyPulseq

```python
seq.add_block(rf, gz, pp.make_delay(TR))

if pp.calc_duration(rf, gz) > TR:
    raise ValueError("events are longer than TR")
seq.add_block(rf, gz, pp.make_delay(TR))
```

:::

`Delay` and `Duration` are construction helpers: they update block `DUR`; they
are not stored as RF, gradient, ADC, or extension events.

## Scanner and Raster Times

Use `Sequence(sys)` when export checks should use a scanner. `write_seq` checks
raster timing and event-local hardware limits from `seq.DEF`, or from `sys` when
passed. Plain `Sequence()` uses Pulseq file rasters and non-limiting hardware
limits.

:::tabs

== KomaMRI

```julia
sys = Scanner(Gmax=40e-3, Smax=150)
seq = Sequence(sys)
@addblock seq += (rf, z=gz)
write_seq(seq, "sequence.seq")  # checks raster and hw limits in seq.DEF
```

== MATLAB Pulseq

```matlab
sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s');
seq = mr.Sequence(sys);
seq.addBlock(rf, gz)
seq.write('sequence.seq')
```

== PyPulseq

```python
system = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s')
seq = pp.Sequence(system=system)
seq.add_block(rf, gz)
seq.write('sequence.seq')
```

:::

## Add Blocks in Loops

Use `@addblocks` around loops that append to a sequence. Inside the macro,
`seq += ...` appends in place.

:::tabs

== KomaMRI

```julia
@addblocks for ky in 1:Ny
    seq += (rf, z=gz)
    seq += (x=gx, y=phase_blip(ky), adc)
end
```

== MATLAB Pulseq

```matlab
for ky = 1:Ny
    seq.addBlock(rf, gz)
    seq.addBlock(adc, gx, phaseBlip(ky))
end
```

== PyPulseq

```python
for ky in range(Ny):
    seq.add_block(rf, gz)
    seq.add_block(adc, gx, phase_blip(ky))
```

:::

Do not use plain `seq += readout` in long loops. In Julia it means
`seq = seq + readout`; Koma's `+` returns a copied sequence so reused events do not share
mutable events. That is safe, but can be more than 800x slower in long loops.

## Reuse Named Sequences

Name a small part of the pulse program as a normal `Sequence`, then append it.
For example, `readout` below is a one-block `Sequence`, not a new event type.

:::tabs

== KomaMRI

```julia
readout = Sequence()
@addblock readout += (ADC(num_readout_samples, adc_duration, ζ), x=Grad(G_readout, T_readout, ζ))

prephaser = Sequence()
@addblock prephaser += (x=Grad(Gx_prephaser, T_prephaser, ζ_prephaser), y=Grad(Gy_prephaser, T_prephaser, ζ_prephaser))

@addblock seq += prephaser + readout
```

== MATLAB Pulseq

```matlab
readout = {gx_readout, adc};
prephaser = {gx_prephaser, gy_prephaser};
addLine(seq, prephaser, readout)

function addLine(seq, prephaser, readout)
    seq.addBlock(prephaser{:})
    seq.addBlock(readout{:})
end
```

== PyPulseq

```python
readout = [gx_readout, adc]
prephaser = [gx_prephaser, gy_prephaser]

def add_line(seq, prephaser, readout):
    seq.add_block(*prephaser)
    seq.add_block(*readout)

add_line(seq, prephaser, readout)
```

:::

Append named sequences in loops with `@addblocks`:

:::tabs

== KomaMRI

```julia
@addblocks for ky in 1:Ny
    seq += rf_preparation + readout(ky)
end
```

== MATLAB Pulseq

```matlab
for ky = 1:Ny
    addRfPreparation(seq)
    addReadout(seq, ky)
end
```

== PyPulseq

```python
for ky in range(Ny):
    add_rf_preparation(seq)
    add_readout(seq, ky)
```

:::

## Transform Sequences

Koma defines arithmetic on `Sequence` values. Each operation returns a copy, so
the original sequence can be reused. Inside `@addblock`, left-multiplying a block
tuple applies the same operation to that one-block sequence.

### Scale Gradients

Use real scalars to scale gradient amplitudes.

:::tabs

== KomaMRI

```julia
@addblock seq += 0.5 * (x=gx, adc)
```

== MATLAB Pulseq

```matlab
seq.addBlock(mr.scaleGrad(gx, 0.5), adc)
```

== PyPulseq

```python
seq.add_block(pp.scale_grad(gx, 0.5), adc)
```

:::

### Rotate Gradients

`real_matrix * sequence` mixes gradient axes and leaves RF and ADC unchanged.

:::tabs

== KomaMRI

```julia
@addblock seq += rotz(π / 6) * (x=gx, adc)
```

== MATLAB Pulseq

```matlab
seq.addBlock(mr.rotate('z', pi/6, gx, adc))
```

== PyPulseq

```python
seq.add_block(*pp.rotate(gx, adc, angle=pi / 6, axis='z'))
```

:::

### Phase RF and ADC

`complex_scalar * sequence` phase-shifts RF and ADC. Gradients are unchanged.
Use `cis(ϕ)` for `exp(im * ϕ)`.

:::tabs

== KomaMRI

```julia
phase = cis(π / 2)  # exp(im * π / 2)
@addblock seq += phase * (rf, z=gz) + phase * (x=gx, adc)
```

== MATLAB Pulseq

```matlab
phase = pi/2;
rf_phase = rf;
adc_phase = adc;
rf_phase.phaseOffset = rf.phaseOffset + phase;
adc_phase.phaseOffset = adc.phaseOffset + phase;

seq.addBlock(rf_phase, gz)
seq.addBlock(gx, adc_phase)
```

== PyPulseq

```python
from copy import copy

phase = pi / 2
rf_phase = copy(rf)
adc_phase = copy(adc)
rf_phase.phase_offset += phase
adc_phase.phase_offset += phase

seq.add_block(rf_phase, gz)
seq.add_block(gx, adc_phase)
```

:::

## Build Radial Spokes

Rotate the readout block for each spoke:

:::tabs

== KomaMRI

```julia
@addblocks for spoke in 0:Nspokes-1
    θ = π * spoke / Nspokes
    seq += (rf, z=gz) + rotz(θ) * (x=gx, adc)
end
```

== MATLAB Pulseq

```matlab
for spoke = 0:Nspokes-1
    theta = pi * spoke / Nspokes;
    seq.addBlock(rf, gz)
    seq.addBlock(mr.rotate('z', theta, gx, adc))
end
```

== PyPulseq

```python
for spoke in range(Nspokes):
    theta = pi * spoke / Nspokes
    seq.add_block(rf, gz)
    seq.add_block(*pp.rotate(gx, adc, angle=theta, axis='z'))
```

:::

## Append Semantics

Use this section when you need to reason about copies or performance.
Tuple terms become one new block through `addblock!`. Existing `Sequence` terms
append all of their blocks through `append!`.

:::tabs

== Source

```julia
@addblock seq += (rf, z=gz) + (x=gx, adc)
```

== Conceptual Calls

```julia
addblock!(seq, rf; z=gz)      # copies rf and gz into a new block, then appends it
addblock!(seq, adc; x=gx)     # copies adc and gx into a new block, then appends it
```

:::

`@addblocks` applies the same rewrite to `+=` expressions in its scope when the
left-hand side is a `Sequence`.

# Create Your Own Sequence

A Koma `Sequence` is made of smaller sequence blocks. A block is a `Sequence` of
length 1: it stores an RF pulse, one gradient per axis, an ADC event, a block
duration, and one or more extensions such as labels, triggers, or rotations.

Use `@addblock` to add one block to a sequence:

:::tabs

== KomaMRI

```julia
seq = Sequence()  # or Sequence(sys) for export checks; see Scanner And Raster Times
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

Everything after `+=` can also compose blocks and reusable chunks. That makes
several Pulseq-style `addBlock` calls one Koma statement:

:::tabs

== KomaMRI

```julia
@addblock seq += (rf, z=gz) + readout(ky) + (Delay(TR), LabelInc(1, "ECO"))
```

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz)
addReadout(seq, ky)
seq.addBlock(mr.makeDelay(TR), mr.makeLabel('INC', 'ECO', 1))
```

== PyPulseq

```python
seq.add_block(rf, gz)
add_readout(seq, ky)
seq.add_block(pp.make_delay(TR), pp.make_label('ECO', 'INC', 1))
```

:::

Unlike MATLAB Pulseq and PyPulseq, a Koma `Grad` does not remember whether it is
an x, y, or z gradient. You choose the axis when adding it to the sequence with
`x=`, `y=`, or `z=`. RF, ADC, and extensions are written normally.

The same block in MATLAB Pulseq or PyPulseq is written by passing the gradient
object for the desired channel directly:

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

Use `Delay(T)` to set a minimum block duration. It is a construction helper that
updates the block `DUR`; it is not stored as an RF, gradient, ADC, or extension.

:::tabs

== KomaMRI

```julia
@addblock seq += (rf, Delay(TR), z=gz)
```

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz, mr.makeDelay(TR))
```

== PyPulseq

```python
seq.add_block(rf, gz, pp.make_delay(TR))
```

:::

Use `Duration(T)` when the block must last exactly `T` seconds. This also sets
`DUR`, but errors if the RF, gradients, or ADC are longer than `T`. MATLAB Pulseq
has the same idea with a numeric first argument to `addBlock`. PyPulseq does not
have a single public `add_block` argument for exact block duration; check the
event duration and add a delay event.

:::tabs

== KomaMRI

```julia
@addblock seq += (rf, Duration(TR), z=gz)
```

== MATLAB Pulseq

```matlab
seq.addBlock(TR, rf, gz)
```

== PyPulseq

```python
if pp.calc_duration(rf, gz) > TR:
    raise ValueError("events are longer than TR")
seq.add_block(rf, gz, pp.make_delay(TR))
```

:::

`@addblock` is a Julia macro: it rewrites the next `seq += ...` expression before
it runs. It does not create sequences; make reusable chunks explicitly, then use
`@addblock` to fill them.

:::tabs

== KomaMRI

```julia
excitation = Sequence()
@addblock excitation += (rf, z=gz)
```

== MATLAB Pulseq

```matlab
excitation = mr.Sequence();
excitation.addBlock(rf, gz)
```

== PyPulseq

```python
excitation = pp.Sequence()
excitation.add_block(rf, gz)
```

:::

## Scanner And Raster Times

Use `Sequence(sys)` when the sequence should carry scanner raster and hardware
metadata. `@addblock` only builds sequence blocks by default. `write_seq` checks
raster timing and event-local hardware limits using `seq.DEF`, or using `sys`
when passed. With plain `Sequence()`, the stored hardware limits are
non-limiting; use `Sequence(sys)` or `write_seq(seq, filename; sys)` for real
scanner checks. Hardware checks currently cover RF amplitude, gradient
amplitude/slew, and ADC dwell; they do not enforce RF ring-down or RF/ADC dead
times.

:::tabs

== KomaMRI

```julia
sys = Scanner(Gmax=40e-3, Smax=150)
seq = Sequence(sys)
@addblock seq += (rf, z=gz)
write_seq(seq, "sequence.seq")
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

`Sequence()` stores Pulseq-style file defaults: `10e-6` seconds for block and
gradient rasters, `1e-6` seconds for RF, and `100e-9` seconds for ADC. `Scanner()`
has scanner/simulation defaults; for example its ADC raster is `2e-6` seconds.
Use `Sequence(sys)` or `write_seq(seq, filename; sys)` when the scanner raster
times should be authoritative.

## Multiple Blocks

For a longer sequence, you will usually add blocks inside a loop. Koma provides
`@addblocks` for this: it lets you write many `seq += ...` lines, and turns them
into efficient block appends.

:::tabs

== KomaMRI

```julia
@addblocks for ky in 1:Ny
    seq += (rf, z=gz)
    seq += (adc, x=gx, y=phase_blip(ky))
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

Use `@addblocks` in loops instead of plain `seq += chunk`. In Julia, `seq += chunk`
is lowered to `seq = seq + chunk`. Koma defines `+` to return a fresh copied
sequence on purpose, so sequence composition is safe and reused chunks do not
share mutable events. In long loops that means the accumulated sequence is copied
again and again; in a 1,000-block loop this can be more than 800x slower.

Inside `@addblocks`, the same rules apply: each tuple is one block and `+`
appends chunks in order.

## Reusable Chunks

You can also give a block or group of blocks a name. This is useful for readouts,
prephasers, refocusing modules, or any piece of sequence you want to reuse.
MATLAB Pulseq and PyPulseq usually do this with normal variables, cell arrays,
lists, or helper functions rather than sequence chunks.

:::tabs

== KomaMRI

```julia
readout = Sequence()
@addblock readout += (ADC(num_readout_samples, adc_duration, ζ), x=Grad(G_readout, T_readout, ζ))

prephaser = Sequence()
@addblock prephaser += (x=Grad(Gx_prephaser, T_prephaser, ζ_prephaser), y=Grad(Gy_prephaser, T_prephaser, ζ_prephaser))

line = prephaser + readout
@addblock seq += line
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

When a chunk is appended, Koma copies its events. That means you can safely reuse
the same chunk many times without later edits to one block changing the others.

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

## Gradient Axes

Named tuples make multi-axis gradients easy to pass around. Here `G_sliceselect`
is one small object with `x`, `y`, and `z` fields.

:::tabs

== KomaMRI

```julia
rf_event = RF(rf_waveform, rf_intervals, Δf, rf_delay)
G_sliceselect = (x=Grad(Gx_sliceselect, T, ζ), y=Grad(Gy_sliceselect, T, ζ), z=Grad(Gz_sliceselect, T, ζ))
G_rewinder = (x=Grad(Gx_rewinder, T_rewinder, ζ), y=Grad(Gy_rewinder, T_rewinder, ζ), z=Grad(Gz_rewinder, T_rewinder, ζ))

excitation = Sequence()
@addblock excitation += (rf_event; G_sliceselect...) + (; G_rewinder...)
```

== MATLAB Pulseq

```matlab
G_sliceselect = {gx_sliceselect, gy_sliceselect, gz_sliceselect};
G_rewinder = {gx_rewinder, gy_rewinder, gz_rewinder};

seq.addBlock(rf, G_sliceselect{:})
seq.addBlock(G_rewinder{:})
```

== PyPulseq

```python
G_sliceselect = [gx_sliceselect, gy_sliceselect, gz_sliceselect]
G_rewinder = [gx_rewinder, gy_rewinder, gz_rewinder]

seq.add_block(rf, *G_sliceselect)
seq.add_block(*G_rewinder)
```

:::

You can also splice variable block arguments.

:::tabs

== KomaMRI

```julia
contents = (Delay(TR), LabelInc(1, "ECO"))
@addblock seq += (contents..., x=gx, y=gy, z=gz)
```

== MATLAB Pulseq

```matlab
contents = {mr.makeDelay(TR), mr.makeLabel('INC', 'ECO', 1)};
seq.addBlock(contents{:}, gx, gy, gz)
```

== PyPulseq

```python
contents = [pp.make_delay(TR), pp.make_label('ECO', 'INC', 1)]
seq.add_block(*contents, gx, gy, gz)
```

:::

## Sequence Arithmetic

Julia's multiple dispatch lets Koma define useful operations on sequence chunks.
Each operation returns a copy, so reusable chunks stay independent.

### Gradient Scaling

`real_scalar * sequence` scales gradients and leaves RF and ADC unchanged.

:::tabs

== KomaMRI

```julia
readout = Sequence()
@addblock readout += (adc, x=gx)

@addblock seq += 0.5 * readout
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

### Gradient Rotation

`real_matrix * sequence` mixes gradient axes and leaves RF and ADC unchanged. A
rotation around z is written with `rotz(ϕ)`:

:::tabs

== KomaMRI

```julia
readout = Sequence()
@addblock readout += (adc, x=gx)

ϕ = π / 6
@addblock seq += rotz(ϕ) * readout
```

== MATLAB Pulseq

```matlab
phi = pi/6;
seq.addBlock(mr.rotate('z', phi, gx, adc))
```

== PyPulseq

```python
phi = pi / 6
seq.add_block(*pp.rotate(gx, adc, angle=phi, axis='z'))
```

:::

### RF And ADC Phase

`complex_scalar * sequence` phase-shifts RF and ADC. Gradients are unchanged.
Use `cis(ϕ)` for a pure phase, equivalent to `exp(im * ϕ)`.

:::tabs

== KomaMRI

```julia
excitation = Sequence()
@addblock excitation += (rf, z=gz)

readout = Sequence()
@addblock readout += (adc, x=gx)

phase = cis(π / 2)  # exp(im * π / 2)
@addblock seq += phase * (excitation + readout)
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

## Radial Readouts

A rotation matrix makes radial readouts compact: define one readout module, then
rotate its gradients for each spoke.

:::tabs

== KomaMRI

```julia
@addblocks for spoke in 0:Nspokes-1
    θ = π * spoke / Nspokes
    seq += excitation + rotz(θ) * readout
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

## What The Macro Does

The macro is syntax for explicit block appends. Conceptually:

:::tabs

== Source

```julia
@addblock seq += (rf, z=gz) + readout(ky)
```

== Equivalent

```julia
addblock!(seq, rf; z=gz)
append!(seq, copy(readout(ky)))
```

:::

Named chunks are normal `Sequence` values. The macro only appends blocks to them:

:::tabs

== Source

```julia
readout = Sequence()
@addblock readout += (adc, x=gx)
@addblock seq += rf_preparation + readout
```

== Equivalent

```julia
readout = Sequence()
addblock!(readout, adc; x=gx)
append!(seq, copy(rf_preparation))
append!(seq, copy(readout))
```

:::

The explicit `copy` in the equivalent code is the important behavior: incoming
RF, gradient, ADC, and extension events are not shared with the destination
sequence.

Non-`Sequence` `+=` expressions keep their normal Julia meaning.

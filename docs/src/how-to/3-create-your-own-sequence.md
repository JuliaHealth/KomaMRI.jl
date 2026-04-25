# Create Your Own Sequence

A Koma `Sequence` is made of smaller sequence blocks. A block is a `Sequence` of
length 1: it stores an RF pulse, one gradient per axis, an ADC event, a block
duration, and extensions such as labels or triggers.

Use `@addblock` to add one block:

:::tabs

== KomaMRI

```julia
seq = Sequence()
@addblock seq += (rf, z=gz)
```

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz)
```

== PyPulseq

```python
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
raster timing and hardware limits before writing.

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

Set Pulseq raster definitions on the sequence only when you need explicit file
rasters without passing a `Scanner` to `write_seq`:

```julia
seq.DEF["BlockDurationRaster"] = sys.DUR_Δt
seq.DEF["GradientRasterTime"] = sys.GR_Δt
seq.DEF["RadiofrequencyRasterTime"] = sys.RF_Δt
seq.DEF["AdcRasterTime"] = sys.ADC_Δt
```

Koma's default `Scanner` matches MATLAB Pulseq for block duration, gradient, and
RF rasters. The ADC raster is different: Koma defaults to `2e-6` seconds, while
MATLAB Pulseq defaults to `100e-9` seconds. Use `write_seq(seq, filename; sys)`
when you want the scanner raster times to be authoritative.

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

Julia's multiple dispatch lets packages define what operators mean for their own
types. Koma uses this for sequence chunks:

- `real_scalar * sequence` copies the sequence and scales gradient amplitudes.
- `real_matrix * sequence` copies the sequence and mixes the gradient axes.
- `complex_scalar * sequence` copies the sequence, phase-shifts RF and ADC, and
  leaves gradients unchanged.

```julia
@addblock seq += 0.5 * readout(ky)  # half gradient amplitude
```

A `3 × 3` matrix mixes the `x`, `y`, and `z` gradient axes. This is useful for
rotating a sequence module:

```julia
θ = π / 6
Rz = [cos(θ) -sin(θ) 0
      sin(θ)  cos(θ) 0
      0       0      1]

@addblock seq += Rz * readout(ky)
```

## Radial Readouts

A rotation matrix makes radial readouts compact: define one readout module, then
rotate its gradients for each spoke.

:::tabs

== KomaMRI

```julia
@addblocks for spoke in 0:Nspokes-1
    θ = π * spoke / Nspokes
    Rz = [cos(θ) -sin(θ) 0
          sin(θ)  cos(θ) 0
          0       0      1]
    seq += Rz * (excitation + readout)
end
```

== MATLAB Pulseq

```matlab
for spoke = 0:Nspokes-1
    theta = pi * spoke / Nspokes;
    seq.addBlock(rf, gz)
    addRotatedReadout(seq, theta)
end
```

== PyPulseq

```python
for spoke in range(Nspokes):
    theta = pi * spoke / Nspokes
    seq.add_block(rf, gz)
    add_rotated_readout(seq, theta)
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

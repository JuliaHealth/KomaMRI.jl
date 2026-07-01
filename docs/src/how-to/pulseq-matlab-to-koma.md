# Pulseq MATLAB to Koma Translation Tables

Compact lookup tables for humans and LLMs translating MATLAB Pulseq code to Koma.
For the block-construction syntax behind `@addblock` and `@addblocks`, see
[Build Sequences with `@addblock`](3-create-your-own-sequence.md).
For a complete translated sequence, see
[Building and Exporting a Pulseq GRE Sequence](../tutorial/gen-09-PulseqGradientEcho.md).

## Conventions

PulseDesigner `make_*` constructor inputs follow Pulseq units for plain numbers.
That is the external sequence-design API. Koma event objects and `Scanner` fields
store physical SI values internally, so accessors such as `area(g)` return SI
even when the constructor input matched Pulseq.

Unitful is per argument, not per call. The same constructor call can mix
Unitful quantities, plain Pulseq numerics, and dimensionless options.

```julia
gx = make_trapezoid(; flat_area=Nx / FOV, flat_time=6.4u"ms", sys)

rf, gz, gzr = make_sinc_pulse(
    90u"deg"; duration=3u"ms", slice_thickness=thickness,
    time_bw_product=4, apodization=0.5, sys,
)
```

## Sequence setup

:::tabs

== MATLAB Pulseq

```matlab
addpath(genpath('path/to/pulseq/matlab'))  % Import Pulseq

sys = mr.opts( ...                         % Define scanner
    'MaxGrad', 32, 'GradUnit', 'mT/m', ...
    'MaxSlew', 130, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 30e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6);

seq = mr.Sequence(sys)                     % Create sequence
seq.setDefinition('FOV', [fov fov dz])     % Set definition
seq.write('seq.seq')                       % Write file
```

== Koma

```julia
using KomaMRI
using KomaMRI.PulseDesigner

sys = Scanner(                         # Define scanner
    Gmax=32e-3, Smax=130.0,
    RF_ring_down_time=30e-6,
    RF_dead_time=100e-6,
    ADC_dead_time=10e-6,
)

seq = Sequence(sys)                    # Create sequence
seq.DEF["FOV"] = [fov, fov, dz]        # Set definition
write_seq(seq, "seq.seq"; sys)         # Write file
```

== Koma + Unitful

```julia
using KomaMRI
using KomaMRI.PulseDesigner
using Unitful

sys = Scanner(                         # Define scanner
    Gmax=32u"mT/m", Smax=130u"T/m/s",
    RF_ring_down_time=30u"μs",
    RF_dead_time=100u"μs",
    ADC_dead_time=10u"μs",
)

seq = Sequence(sys)                         # Create sequence
seq.DEF["FOV"] = [FOV, FOV, dz] .|> to_SI   # Set definition
write_seq(seq, "seq.seq"; sys)              # Write file
```

:::

## Blocks and axes

:::tabs

== MATLAB Pulseq

```matlab
seq.addBlock(rf, gz)
seq.addBlock(gx, adc)
seq.addBlock(gxPre, gyPre, gzReph)
seq.addBlock(mr.makeLabel('SET', 'LIN', ky))
seq.addBlock(gxPre, mr.scaleGrad(gyPre, peScale))
if hasAdc
    seq.addBlock(gx, adc)
else
    seq.addBlock(gx)
end
```

== Koma @addblock

```julia
@addblock seq += (rf, z=gz)
@addblock seq += (x=gx, adc)
@addblock seq += (x=gx_pre, y=gy_pre, z=gz_reph)
@addblock seq += make_label(:SET, :LIN, ky)
@addblock seq += (x=gx_pre, y=pe_scale * gy_pre)
@addblock seq += (x=gx, has_adc ? adc : nothing)
```

== Koma @addblocks

```julia
@addblocks begin
    seq += (rf, z=gz)
    seq += (x=gx, adc)
    seq += (x=gx_pre, y=gy_pre, z=gz_reph)
    seq += make_label(:SET, :LIN, ky)
    seq += (x=gx_pre, y=pe_scale * gy_pre)
    seq += (x=gx, has_adc ? adc : nothing)
end
```

:::

## Timing helpers

Use `sys` when translating MATLAB formulas that depend on sample-edge,
ring-down, dead-time, or block duration behavior.

:::tabs

== MATLAB Pulseq

```matlab
ceil(t / raster) * raster             % Round up to raster
round(t / raster) * raster            % Round to raster

mr.calcDuration(rf)                   % RF duration
mr.calcDuration(adc)                  % ADC duration
mr.calcDuration(rf, gz)               % Block duration
mr.calcRfCenter(rf)                   % RF center
adc.dwell                             % ADC dwell
```

== Koma

```julia
ceil_to_raster(t, sys.GR_Δt)          # Round up to raster
round_to_raster(t, sys.GR_Δt)         # Round to raster

dur(rf, sys)                          # RF duration with sample-edge timing and ring-down
dur(adc, sys)                         # ADC duration with dwell-edge timing and dead time
dur(seq[1], sys)                      # One-block duration with sys-aware timing
rf_center(rf, sys)                    # RF center with sample-edge timing
dwell(adc, sys)                       # ADC dwell

delay(rf, sys)                        # RF delay with sample-edge timing
delay(adc, sys)                       # ADC delay with dwell-edge timing
```

:::

## Trapezoidal gradients

Plain PulseDesigner gradient numerics use Pulseq units. Koma event accessors
such as `area(gx)` return physical SI values.

:::tabs

== MATLAB Pulseq

```matlab
fov = 0.22;
Nx = 64;
sliceThickness = 3e-3;
deltaK = 1/fov;
readoutTime = 6.4e-3;
preTime = 2e-3;
blipDuration = 0.5e-3;
testAmplitude = 10e3;
testRiseTime = 1e-3;
testFlatTime = 10e-3;

gx = mr.makeTrapezoid('x', 'FlatArea', Nx*deltaK, 'FlatTime', readoutTime, 'system', sys)
gxEpi = mr.makeTrapezoid('x', 'Amplitude', Nx*deltaK/readoutTime, 'FlatTime', readoutTime, 'system', sys)
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', preTime, 'system', sys)
gyBlip = mr.makeTrapezoid('y', 'Area', deltaK, 'Duration', blipDuration, 'system', sys)
gzSpoil = mr.makeTrapezoid('z', 'Area', 4/sliceThickness, 'system', sys)
gTest = mr.makeTrapezoid('x', 'Amplitude', testAmplitude, 'riseTime', testRiseTime, 'flatTime', testFlatTime, 'fallTime', 2*testRiseTime, 'system', sys)
```

== Koma

```julia
fov = 0.22
Nx = 64
slice_thickness = 3e-3
delta_k = 1 / fov
readout_time = 6.4e-3
pre_time = 2e-3
blip_duration = 0.5e-3
test_amplitude = 10e3
test_rise_time = 1e-3
test_flat_time = 10e-3

gx = make_trapezoid(; flat_area=Nx * delta_k, flat_time=readout_time, sys)
gx_epi = make_trapezoid(; amplitude=Nx * delta_k / readout_time, flat_time=readout_time, sys)
gx_pre = make_trapezoid(; area=-area(gx) * γ / 2, duration=pre_time, sys)
gy_blip = make_trapezoid(; area=delta_k, duration=blip_duration, sys)
gz_spoil = make_trapezoid(; area=4 / slice_thickness, sys)
g_test = make_trapezoid(;
    amplitude=test_amplitude, rise_time=test_rise_time,
    flat_time=test_flat_time, fall_time=2test_rise_time, sys,
)
```

== Koma + Unitful

```julia
FOV = 220u"mm"
Nx = 64
slice_thickness = 3u"mm"
delta_k = 1 / FOV
readout_time = 6.4u"ms"
pre_time = 2u"ms"
blip_duration = 0.5u"ms"
test_amplitude = 10_000u"Hz/m"
test_rise_time = 1u"ms"
test_flat_time = 10u"ms"

gx = make_trapezoid(; flat_area=Nx * delta_k, flat_time=readout_time, sys)
gx_epi = make_trapezoid(; amplitude=Nx * delta_k / readout_time, flat_time=readout_time, sys)
gx_pre = make_trapezoid(; area=-area(gx) / 2 * u"T*s/m", duration=pre_time, sys)
gy_blip = make_trapezoid(; area=delta_k, duration=blip_duration, sys)
gz_spoil = make_trapezoid(; area=4 / slice_thickness, sys)
g_test = make_trapezoid(;
    amplitude=test_amplitude, rise_time=test_rise_time,
    flat_time=test_flat_time, fall_time=2test_rise_time, sys,
)
```

:::

## Extended gradients and gradient sums

:::tabs

== MATLAB Pulseq

```matlab
gyBlipUp = mr.makeExtendedTrapezoid('y', 'times', timesUp, 'amplitudes', ampsUp, 'system', sys)
gz180n = mr.makeExtendedTrapezoid('z', 'times', gz180Times, 'amplitudes', gz180Amps, 'system', sys)
[gzr1, times, amplitudes] = mr.makeExtendedTrapezoidArea('z', g0, g1, A, sys)
gSum = mr.addGradients({g1, g2}, sys)
```

== Koma

```julia
gy_blipup = make_extended_trapezoid(times_up, amps_up; sys)
gz180n = make_extended_trapezoid(gz180_times, gz180_amps; sys)
gzr1 = make_extended_trapezoid_area(g0, g1, A; sys)
g_sum = g1 + g2
```

:::

## RF pulses

:::tabs

== MATLAB Pulseq

```matlab
rf = mr.makeBlockPulse(pi/2, 'Duration', 300e-6, 'system', sys)

[rf, gz, gzr] = mr.makeSincPulse(pi/2, 'Duration', 3e-3, 'SliceThickness', thickness, 'system', sys)
[rf, gz, gzr] = mr.makeSincPulse(pi/2, 'Duration', 3e-3, 'SliceThickness', thickness, 'apodization', 0.5, 'timeBwProduct', 4, 'system', sys)
rfInv = mr.makeSincPulse(pi, 'Duration', 10e-3, 'timeBwProduct', 8, 'system', sys)

[rfFat, ~, ~] = mr.makeGaussPulse(pi/2, 'Duration', 8e-3, 'Bandwidth', abs(fatOffset), 'freqOffset', fatOffset, 'system', sys)
[rf, ~, ~] = mr.makeArbitraryRf(rfSignal, flipAngle, 'Dwell', dwell, 'Delay', rfDelay, 'system', sys)
```

== Koma

```julia
rf = make_block_pulse(π / 2; duration=300e-6, sys)

rf, gz, gzr = make_sinc_pulse(π / 2; duration=3e-3, slice_thickness=thickness, sys)
rf, gz, gzr = make_sinc_pulse(π / 2; duration=3e-3, slice_thickness=thickness, apodization=0.5, time_bw_product=4, sys)
rf_inv, _, _ = make_sinc_pulse(π; duration=10e-3, time_bw_product=8, sys)

rf_fat, _, _ = make_gauss_pulse(π / 2; duration=8e-3, bandwidth=abs(fat_offset), freq_offset=fat_offset, sys)
rf, _, _ = make_arbitrary_rf(rf_signal, flip_angle; dwell, delay=rf_delay, sys)
```

== Koma + Unitful

```julia
rf = make_block_pulse(90u"deg"; duration=300u"μs", sys)

rf, gz, gzr = make_sinc_pulse(90u"deg"; duration=3u"ms", slice_thickness=thickness, sys)
rf, gz, gzr = make_sinc_pulse(90u"deg"; duration=3u"ms", slice_thickness=thickness, apodization=0.5, time_bw_product=4, sys)
rf_inv, _, _ = make_sinc_pulse(180u"deg"; duration=10u"ms", time_bw_product=8, sys)

rf_fat, _, _ = make_gauss_pulse(90u"deg"; duration=8u"ms", bandwidth=abs(fat_offset), freq_offset=fat_offset, sys)
rf, _, _ = make_arbitrary_rf(rf_signal, flip_angle; dwell, delay=rf_delay, sys)
```

:::

## ADC, delays, labels, and extensions

:::tabs

== MATLAB Pulseq

```matlab
readoutTime = 6.4e-3;
adcDelay = gx.riseTime;
adcDwell = 4e-6;

adc = mr.makeAdc(Nx, 'Duration', readoutTime, 'Delay', adcDelay, 'system', sys)
adc = mr.makeAdc(adcSamples, 'Dwell', adcDwell, 'Delay', adcDelay, 'system', sys)

delayTE = mr.makeDelay(20e-3)
label = mr.makeLabel('SET', 'LIN', ky)
trig = mr.makeTrigger('physio1', 'duration', 2000e-6)
osc = mr.makeDigitalOutputPulse('osc0', 'duration', 100e-6)

rot = mr.makeRotation(phi)
seq.addBlock(mr.rotate('z', phi, gx), adc)
```

== Koma

```julia
readout_time = 6.4e-3
adc_delay = gx.rise
adc_dwell = 4e-6

adc = make_adc(Nx; duration=readout_time, delay=adc_delay, sys)
adc = make_adc(adc_samples; dwell=adc_dwell, delay=adc_delay, sys)

delay_te = make_delay(20e-3)
label = make_label(:SET, :LIN, ky)
trig = make_trigger(:physio1; duration=2000e-6, sys)
osc = make_digital_output_pulse(:osc0; duration=100e-6, sys)

rot = make_rotation(phi)
@addblock seq += rotz(phi) * (x=gx, adc)
```

== Koma + Unitful

```julia
readout_time = 6.4u"ms"
adc_delay = gx.rise * u"s"
adc_dwell = 4u"μs"

adc = make_adc(Nx; duration=readout_time, delay=adc_delay, sys)
adc = make_adc(adc_samples; dwell=adc_dwell, delay=adc_delay, sys)

delay_te = make_delay(20u"ms")
label = make_label(:SET, :LIN, ky)
trig = make_trigger(:physio1; duration=2000u"μs", sys)
osc = make_digital_output_pulse(:osc0; duration=100u"μs", sys)

rot = make_rotation(phi)
@addblock seq += rotz(phi) * (x=gx, adc)
```

:::

## Not-yet-mapped Pulseq features

These Pulseq features are not direct PulseDesigner APIs or full parity paths
yet. Keep them explicit in translations instead of hiding them behind one-off
local wrappers.

| Pulseq feature | Koma status |
|---|---|
| `mr.makeRfShim(...)` | RF shimming extension is not implemented in Pulseq read/write. |
| `mr.makeSoftDelay(...)` | Soft-delay extension is not implemented in Pulseq read/write. |
| `mr.makeSLRpulse(...)` | No PulseDesigner constructor; port the waveform explicitly if needed. |
| `mr.addRamps(...)` | No ramp-optimization helper; port the final gradient waveform explicitly. |
| `mr.traj2grad(...)` | No direct helper; derive the gradient waveform explicitly. |
| ADC `phaseModulation` | Koma `ADC` stores one phase offset; ADC phase shapes are not strict parity yet. |
| `mr.align(...)` | No direct clone helper; use explicit delays or block `Duration(...)`. |
| `mr.splitGradientAt(...)` | No public clone helper; rebuild the required pieces explicitly when needed. |

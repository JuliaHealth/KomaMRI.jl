# Pulseq Event Mapping

KomaMRI uses two representations:

- [`Sequence`](4-sequence.md): the runtime representation used by simulation and plotting.
  It stores repeated block events as [`RF`](5-seq-events.md#rf),
  [`Grad`](5-seq-events.md#gradient), and [`ADC`](5-seq-events.md#adc) objects.
  See [Sequence](4-sequence.md) for the sequence model.
- [`PulseqSequenceData`](@ref): the Pulseq-native intermediate representation
  used by the reader and writer. It stores `[BLOCKS]` rows, event libraries,
  definitions, the Pulseq version, and the optional signature.

API references: [`read_seq`](@ref), [`read_seq_data`](@ref), [`write_seq`](@ref),
[`write_seq_data`](@ref), and [`PulseqSequenceData`](@ref).

## Reader

[`read_seq`](@ref) first reads the file into [`PulseqSequenceData`](@ref), then
materializes a Koma `Sequence`. Use [`read_seq_data`](@ref) when you want the
Pulseq block/event-library representation without expanding every block into
Koma events.

```@raw html
<div style="overflow-x:auto; margin: 1rem 0;">
<svg viewBox="0 0 900 220" role="img" aria-label="Pulseq reader pipeline" style="min-width: 720px; width: 100%; height: auto;">
  <defs>
    <marker id="reader-arrow" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto">
      <path d="M0,0 L0,6 L9,3 z" fill="#64748b"></path>
    </marker>
  </defs>
  <rect x="35" y="60" width="135" height="70" rx="12" fill="#eef6ff" stroke="#60a5fa" stroke-width="2"></rect>
  <text x="102" y="90" text-anchor="middle" font-size="16" font-weight="700">.seq file</text>
  <text x="102" y="112" text-anchor="middle" font-size="12">Pulseq sections</text>

  <line x1="170" y1="95" x2="225" y2="95" stroke="#64748b" stroke-width="2" marker-end="url(#reader-arrow)"></line>
  <rect x="235" y="48" width="205" height="94" rx="12" fill="#f8fafc" stroke="#94a3b8" stroke-width="2"></rect>
  <text x="337" y="78" text-anchor="middle" font-size="16" font-weight="700">read_seq_data</text>
  <text x="337" y="101" text-anchor="middle" font-size="12">parse sections + scales</text>
  <text x="337" y="120" text-anchor="middle" font-size="12">legacy fixups + optional signature</text>

  <line x1="440" y1="95" x2="495" y2="95" stroke="#64748b" stroke-width="2" marker-end="url(#reader-arrow)"></line>
  <rect x="505" y="60" width="180" height="70" rx="12" fill="#fff7ed" stroke="#fb923c" stroke-width="2"></rect>
  <text x="595" y="90" text-anchor="middle" font-size="16" font-weight="700">PulseqSequenceData</text>
  <text x="595" y="112" text-anchor="middle" font-size="12">ids + event libraries</text>

  <line x1="685" y1="95" x2="740" y2="95" stroke="#64748b" stroke-width="2" marker-end="url(#reader-arrow)"></line>
  <rect x="750" y="48" width="115" height="94" rx="12" fill="#ecfdf5" stroke="#34d399" stroke-width="2"></rect>
  <text x="807" y="79" text-anchor="middle" font-size="16" font-weight="700">Sequence</text>
  <text x="807" y="101" text-anchor="middle" font-size="12">decode once</text>
  <text x="807" y="120" text-anchor="middle" font-size="12">expand blocks</text>

  <text x="450" y="185" text-anchor="middle" font-size="13" fill="#64748b">read_seq_data stops at the Pulseq-native representation; read_seq continues through Sequence(data).</text>
</svg>
</div>
```

The reader has three stages:

1. `parse_pulseq_file` reads sections into `PulseqParsedFile`. It applies
   version-specific scales, for example `γ` conversion for RF/gradient
   amplitudes and `μs`/`ns` conversion for timing fields.
2. `pulseq_sequence_data` normalizes legacy details and returns
   [`PulseqSequenceData`](@ref). This is still Pulseq-native: blocks contain
   integer event ids, and event payloads live in libraries.
3. `Sequence(data)` decodes each library entry once, then expands block ids into
   a Koma `Sequence`.

### RF

| Pulseq RF timing | Koma fields | Koma representation |
|---|---|---|
| no RF event id | `RF(0.0, 0.0)` | [`BlockPulseRF`](5-seq-events.md#rf) |
| `time_shape_id == 0` | `A::Vector`, `T::Number`; `delay = pulseq_delay + Δt_rf/2` | [`UniformlySampledRF`](5-seq-events.md#rf-uniformly-sampled-waveform) |
| `time_shape_id == -1` | `A::Vector`, `T::Number`; half-raster compact timing; `delay = pulseq_delay + Δt_rf/2` | [`UniformlySampledRF`](5-seq-events.md#rf-uniformly-sampled-waveform) |
| `time_shape_id > 0` | `A::Vector`, `T::Vector = diff(time_shape) * Δt_rf`; `delay = pulseq_delay + first(time_shape) * Δt_rf` | [`TimeShapedRF`](5-seq-events.md#rf-time-shaped-waveform) |

Pulseq RF events store magnitude and phase as shapes. Therefore, even a constant
RF pulse in a Pulseq file is read as a sampled RF waveform, not as a scalar
`BlockPulseRF`. `BlockPulseRF` appears for empty RF placeholders and Koma-native
scalar RF construction.

Compact RF timing uses Pulseq's center-sampled convention: Pulseq stores
`pulseq_delay` before the implicit first-sample offset, while Koma stores
`delay` to the first RF sample. For compact RF timing, these differ by
`Δt_rf/2`.

### Gradients

| Pulseq gradient form | Koma fields | Koma representation |
|---|---|---|
| no gradient event id | `Grad(0.0, 0.0)` | [`TrapezoidalGrad`](5-seq-events.md#gradient-trapezoidal-waveform) |
| `[TRAP]` event | `A::Number`, `T::Number`, scalar rise/fall | [`TrapezoidalGrad`](5-seq-events.md#gradient-trapezoidal-waveform) |
| `[GRADIENTS]`, `time_shape_id == 0` | `A::Vector`, `T::Number`, `rise = fall = Δt_gr/2` | [`UniformlySampledGrad`](5-seq-events.md#gradient-uniformly-sampled-waveform) |
| `[GRADIENTS]`, `time_shape_id == -1` | `A::Vector`, `T::Number`, half-raster compact timing, `rise = fall = Δt_gr/2` | [`UniformlySampledGrad`](5-seq-events.md#gradient-uniformly-sampled-waveform) |
| `[GRADIENTS]`, `time_shape_id > 0` | `A::Vector`, `T::Vector = diff(time_shape) * Δt_gr`; `delay += first(time_shape) * Δt_gr`; `rise = fall = 0` | [`TimeShapedGrad`](5-seq-events.md#gradient-time-shaped-waveform) |

For Pulseq versions before 1.5, the reader also reconstructs missing arbitrary
gradient `first`/`last` values from neighboring blocks before materializing the
Koma sequence.

### ADC

[`ADC`](5-seq-events.md#adc) events do not have shape-based representations. The timing
convention is the only translation: Koma stores `delay` to the first acquired
sample, while Pulseq stores `pulseq_delay` to the dwell interval edge.

| Pulseq ADC timing | Koma fields |
|---|---|
| no ADC event id | `ADC(0, 0.0)` |
| `[ADC]` event | `N = num`, `T = dwell * (num - 1)`, `delay = pulseq_delay + dwell/2` |

On read, `T` is the time from the first sample to the last sample. A single
sample therefore reads as `T = 0`.

## Writer

[`write_seq`](@ref) converts a Koma `Sequence` into [`PulseqSequenceData`](@ref)
and then emits Pulseq sections. Use [`write_seq_data`](@ref) to get that
Pulseq-native representation without writing a file. `write_seq(data::PulseqSequenceData, filename)`
writes an already prepared representation directly.

```@raw html
<div style="overflow-x:auto; margin: 1rem 0;">
<svg viewBox="0 0 900 220" role="img" aria-label="Pulseq writer pipeline" style="min-width: 720px; width: 100%; height: auto;">
  <defs>
    <marker id="writer-arrow" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto">
      <path d="M0,0 L0,6 L9,3 z" fill="#64748b"></path>
    </marker>
  </defs>
  <rect x="35" y="60" width="135" height="70" rx="12" fill="#ecfdf5" stroke="#34d399" stroke-width="2"></rect>
  <text x="102" y="90" text-anchor="middle" font-size="16" font-weight="700">Sequence</text>
  <text x="102" y="112" text-anchor="middle" font-size="12">runtime events</text>

  <line x1="170" y1="95" x2="225" y2="95" stroke="#64748b" stroke-width="2" marker-end="url(#writer-arrow)"></line>
  <rect x="235" y="42" width="210" height="106" rx="12" fill="#f8fafc" stroke="#94a3b8" stroke-width="2"></rect>
  <text x="340" y="70" text-anchor="middle" font-size="16" font-weight="700">write_seq_data</text>
  <text x="340" y="94" text-anchor="middle" font-size="12">hardware checks</text>
  <text x="340" y="113" text-anchor="middle" font-size="12">raster + quantize</text>
  <text x="340" y="132" text-anchor="middle" font-size="12">timing check + deduplicate</text>

  <line x1="445" y1="95" x2="500" y2="95" stroke="#64748b" stroke-width="2" marker-end="url(#writer-arrow)"></line>
  <rect x="510" y="60" width="180" height="70" rx="12" fill="#fff7ed" stroke="#fb923c" stroke-width="2"></rect>
  <text x="600" y="90" text-anchor="middle" font-size="16" font-weight="700">PulseqSequenceData</text>
  <text x="600" y="112" text-anchor="middle" font-size="12">ids + event libraries</text>

  <line x1="690" y1="95" x2="745" y2="95" stroke="#64748b" stroke-width="2" marker-end="url(#writer-arrow)"></line>
  <rect x="755" y="60" width="110" height="70" rx="12" fill="#eef6ff" stroke="#60a5fa" stroke-width="2"></rect>
  <text x="810" y="90" text-anchor="middle" font-size="16" font-weight="700">.seq file</text>
  <text x="810" y="112" text-anchor="middle" font-size="12">sections + signature</text>

  <text x="450" y="185" text-anchor="middle" font-size="13" fill="#64748b">write_seq(seq; sys) uses sys for checks and rasters; write_seq(data, filename) skips Koma rasterization.</text>
</svg>
</div>
```

The writer has four stages:

1. By default, `check_hw_limits=true` checks hardware limits before writing.
   Without `sys`, this uses metadata in `seq.DEF`; with `sys`, it uses `sys`.
2. `prepare_pulseq_write` copies and quantizes RF, gradient, ADC, and block
   duration timings to the Pulseq rasters. Without `sys`, raster times come from
   `seq.DEF`; with `sys`, the scanner rasters are used.
3. By default, `check_timing=true` checks the quantized copy against the Pulseq
   rasters.
4. `collect_pulseq_assets` and `emit_pulseq` write `[DEFINITIONS]`, `[BLOCKS]`,
   event libraries, shape libraries, extensions, and the signature.

Hardware-limit definitions such as `MaxGrad`, `MaxSlew`, `MaxB1`, and dead times
are Koma check metadata. They are not emitted to the Pulseq `[DEFINITIONS]`
section; the writer emits Pulseq raster definitions and non-internal user
definitions.

### RF

| Koma RF condition | Pulseq output |
|---|---|
| `!is_on(rf)` | no RF event id |
| [`BlockPulseRF`](5-seq-events.md#rf) | `[RF]` event with constant magnitude and phase shapes |
| [`UniformlySampledRF`](5-seq-events.md#rf-uniformly-sampled-waveform) with sample times `0, Δt_rf, 2Δt_rf, ...` and compatible delay | `[RF]` with `time_shape_id = 0`; `pulseq_delay = rf.delay - Δt_rf/2` |
| [`UniformlySampledRF`](5-seq-events.md#rf-uniformly-sampled-waveform) with half-raster sample times and odd sample count | `[RF]` with `time_shape_id = -1`; `pulseq_delay = rf.delay - Δt_rf/2` |
| [`TimeShapedRF`](5-seq-events.md#rf-time-shaped-waveform) | `[RF]` with an explicit time shape |
| otherwise | `[RF]` with an explicit time shape |

`TimeShapedRF` is always written with an explicit time shape. This preserves
vector timing and avoids replacing it with Pulseq's compact RF convention, where
the first RF sample is implicitly offset by `Δt_rf/2`.

Frequency-modulated RF waveforms are not written to Pulseq: `FrequencyModulatedRF`
throws an error in the writer.

### Gradients

| Koma gradient condition | Pulseq output |
|---|---|
| `!is_on(gr)` and zero duration | no gradient event id |
| [`TrapezoidalGrad`](5-seq-events.md#gradient-trapezoidal-waveform) with `first == last == 0` | `[TRAP]` |
| [`TrapezoidalGrad`](5-seq-events.md#gradient-trapezoidal-waveform) with nonzero edge continuity | converted to an edge-timed arbitrary gradient, then `[GRADIENTS]` |
| [`UniformlySampledGrad`](5-seq-events.md#gradient-uniformly-sampled-waveform) with two equal samples and `first == last == 0` | `[TRAP]` |
| [`TimeShapedGrad`](5-seq-events.md#gradient-time-shaped-waveform) with two equal samples, one interval, and `first == last == 0` | `[TRAP]` |
| arbitrary gradient sample times match `(i - 1/2)Δt_gr` | `[GRADIENTS]` with `time_shape_id = 0` |
| arbitrary gradient sample times match `i * Δt_gr/2` with odd sample count | `[GRADIENTS]` with `time_shape_id = -1` |
| otherwise | `[GRADIENTS]` with an explicit time shape |

Unlike RF, a `TimeShapedGrad` may still be compacted if its sample times exactly
match a Pulseq compact gradient timing convention. Gradient compact timing does
not carry the same RF first-sample convention, so this is a representation
optimization.

### ADC

| Koma ADC condition | Pulseq output |
|---|---|
| `!is_on(adc)` | no ADC event id |
| `N == 1` | `[ADC]` with `dwell = T`, `pulseq_delay = delay - dwell/2` |
| `N > 1` | `[ADC]` with `dwell = T / (N - 1)`, `pulseq_delay = delay - dwell/2` |

Pulseq ADC timing follows an edge-versus-sample convention: the Pulseq delay
points to the dwell interval edge, and the first sample is acquired at
`dwell/2`. Koma stores `delay` at the first acquired sample, so the writer
subtracts `dwell/2` and the reader adds it back.

For `N == 1`, the writer uses `T` as the dwell time because there is no
inter-sample interval from which to infer it.

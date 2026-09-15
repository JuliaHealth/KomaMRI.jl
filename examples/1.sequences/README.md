# UI sequence examples

These files use Pulseq 1.5.1, with explicit `Nx`, `Ny`, `Nz`, FOV in metres,
RF uses, and acquisition labels. `Nz=1` describes a 2D encoding matrix;
separate slices are identified by `SLC`, not by RF frequency offsets.

| File | Encoding matrix | ADC readouts × samples | Labels |
| --- | --- | --- | --- |
| `epi_JEMRIS.seq` | 64 × 64 | 64 × 64 | `LIN`, `SLC=0`, `REV` |
| `epi_multislice.seq` | 100 × 100 | 3 × 100 × 100 | `LIN`, `SLC=0:2`, `REV` |
| `epi_ramp.seq` | 64 × 64 | 56 × 84 | `LIN=55:-1:0`, `SLC=0`, `REV` |
| `epi_ramp_fatsat.seq` | 64 × 64 | 56 × 84 | `LIN=55:-1:0`, `SLC=0`, `REV` |
| `epi_se.seq` | 65 × 64 | 64 × 65 | `LIN`, `SLC=0`, `REV` |
| `epi_se_128.seq` | 128 × 128 | 128 × 128 | `LIN`, `SLC=0`, `REV` |
| `ge.seq` | 101 × 100 | 100 × 101 | `LIN`, `SLC=0` |
| `gre_JEMRIS.seq` | 32 × 32 | 32 × 32 | `LIN`, `SLC=0` |
| `slice_profile_z.seq` | 1 × 1 × 256 | 1 × 256 | `SLC=0` |
| `spiral.seq` | 64 × 64 | 1 × 12000 | `SEG=0`, `SLC=0` |

All ADCs are imaging data. `LIN` identifies Cartesian phase-encoding lines;
`REV` marks reversed readouts. The spiral has one continuous ADC event.

`epi_se_128.seq` is single-shot spin-echo EPI with a 256 mm FOV (2 mm in-plane
pixels), 3 mm slice thickness, and approximately 140 ms TE. It uses a 90° sinc
excitation, a 1.2 ms 180° refocusing pulse with symmetric crushers, and 128
alternating readouts. Gradient amplitude, slew rate, and RF amplitude stay within
32 mT/m, 130 T/m/s, and 10 μT respectively.

`slice_profile_z.seq` measures a 1D slice profile: 30° windowed-sinc excitation
(4 ms, time-bandwidth product 8, 5 mm slice), slice rephasing, and frequency
encoding all along z. The readout covers a 40 mm FOV with 256 samples. Use a
uniform, on-resonance phantom spanning that FOV; `Nz=256` is the z readout
matrix, not a partition or slice count. No `PAR` label is needed.

The ramp-sampled EPI examples acquire 56 of 64 phase-encoding lines and 84
samples per ramp-sampled readout. Neither ADC count is the image matrix.
The original odd readout matrices in `epi_se.seq` and `ge.seq` are retained;
Koma currently pads their reconstruction sizes to 66 and 102 respectively.
The fat-saturation pulse is marked as saturation, not another slice excitation.

The upgrade preserves the saved RF/gradient waveforms, block durations, and ADC
sampling times. ADC trajectories are unchanged. Legacy millimetre FOV definitions were
converted to metres; a zero through-plane FOV remains unspecified, not a physical
slice thickness. Three legacy files declare a 1 ns ADC raster to retain their
original non-100-ns-multiple dwell times without rounding.

`Pulseq_create_seq/` contains historical MATLAB generators whose parameters do
not all match these saved acquisitions. The migration annotates the saved
waveforms; rerunning those generators is not an equivalent regeneration.

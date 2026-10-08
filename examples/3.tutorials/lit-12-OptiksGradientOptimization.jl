# # Optimizing a Spiral Readout with OptiksMRI.jl
#
# The time-optimal gradient waveform for a k-space trajectory runs the hardware at its
# limits, and its spectrum can land on the mechanical resonances of the gradient coil.
# [OptiksMRI.jl](https://code.stanford.edu/mcnab-lab/optiksmri.jl) is a pure-Julia port of
# [OPTIKS](https://github.com/mamccready/optiks): it keeps the k-space path fixed and
# optimizes the *timing* along it against a loss you compose from terms such as readout
# duration, slew rate, forbidden frequency bands, and PNS.
#
# In this tutorial we design a variable-density spiral interleaf with OptiksMRI, push its
# spectrum out of a few forbidden frequency bands, and then turn the result into a
# KomaMRI `Sequence` that can be plotted, simulated, or exported to Pulseq.
#
# If you use OPTIKS, please cite: M. A. McCready, X. Cao, K. Setsompop, J. M. Pauly and
# A. B. Kerr, "OPTIKS: Optimized Gradient Properties Through Timing in K-Space," *IEEE
# Transactions on Medical Imaging*, doi:
# [10.1109/TMI.2025.3639398](https://doi.org/10.1109/TMI.2025.3639398).

using KomaMRI
using OptiksMRI

# Both packages export a `plot_kspace` function, so below we call them as
# `OptiksMRI.plot_kspace` and `KomaMRI.plot_kspace`.
#
# ## Scanner and trajectory
#
# We use Koma's default hardware limits: 60 mT/m, 500 T/m/s, and a 10 μs gradient raster.
# OptiksMRI works in SI units, as Koma does, so the limits can be passed straight through.
scanner = Scanner()
lim = scanner.limits;

# The k-space path is one interleaf of a variable-density spiral (`α = 2`) with a 220 mm
# field of view, 1 mm resolution, and an in-plane acceleration of `R = 4`.
# `generic_trajectory` samples it as a matrix with one row per point and one column per
# axis, in units of 1/m.
fov = 220e-3
res = 1e-3
spiral = generic_trajectory(InterleavedSpiral(fov, res; α=2, R=4));

# ## Setting up the optimization
#
# `HardwareOpts` holds the boundary conditions and limits of the waveform: here it starts
# and ends at zero gradient amplitude.
hardware_opts = HardwareOpts(;
    g0   = 0,         # initial gradient amplitude
    gfin = 0,         # final gradient amplitude
    gmax = lim.Gmax,  # gradient amplitude limit
    smax = lim.Smax,  # slew rate limit
    dt   = lim.GR_Δt, # gradient raster
);

# `DesignOpts` composes the loss. We ask for a short readout that respects the slew limit
# while keeping as little energy as possible in the 400–600 Hz and 750–950 Hz bands, and
# above 1750 Hz.
design_opts = DesignOpts(;
    terms=(
        TimeMin(; weight=1e-3),
        SlewLimit(; weight=1e-3),
        FrequencyMin(; weight=1e0, bands=[(400, 600), (750, 950), (1750, Inf)]),
    ),
);

# `SolverOpts` controls the gradient descent. The starting point is the time-optimal
# solution slowed down by `derate`, and `oversample` sets how many optimization variables
# there are per gradient raster step. The AdamW learning rate `lr` usually needs some
# tuning for a new trajectory: here it was adjusted until enough energy was removed from
# the forbidden bands. Set `verbose = true` to watch each loss term during the descent.
solver_opts = SolverOpts(; oversample=8, derate=0.95, lr=2e-3, verbose=false);

# ## Running OptiksMRI
result = optiks(spiral; hardware=hardware_opts, design=design_opts, solver=solver_opts)
result.waveform

# `result.waveform_init` holds the time-optimal initial waveform, so we can compare the
# spectra before and after the optimization. The shaded regions are the forbidden bands.
p1 = plot_spectrum(result.waveform_init, design_opts.terms[3])
#jl display(p1);

#-
p2 = plot_spectrum(result.waveform, design_opts.terms[3])
#jl display(p2);

# The optimized waveform stays within the gradient and slew limits,
p3 = plot_gradient(result.waveform, hardware_opts)
#jl display(p3);

#-
p4 = plot_slew(result.waveform, hardware_opts)
#jl display(p4);

# and it still traces the spiral we asked for.
p5 = OptiksMRI.plot_kspace(result.waveform)
#jl display(p5);

# ## Building a Koma sequence
#
# `result.waveform.T` is the continuous optimal duration, which is generally not a raster
# multiple. The samples in `result.waveform.g` are spaced `hardware_opts.dt` apart, so we
# take the gradient duration from the sample count to land them exactly on the gradient
# raster.
T_grad = (size(result.waveform.g, 1) - 1) * hardware_opts.dt;

# Pulseq opens the ADC window `ADC_dead_time` after the block starts, so we delay the
# gradients by the same amount to keep gradient sample `k` aligned with ADC sample `k`.
D = lim.ADC_dead_time;

# We sample once per gradient raster step. The dwell time must be a multiple of `ADC_Δt`
# (here 10 μs is 5 × 2 μs), and the number of samples a multiple of 4. Koma's `adc.T`
# spans the first to the last sample *centre*, and `adc.delay` is the first sample centre,
# which Pulseq places half a dwell into the window.
dwell = lim.GR_Δt
Nadc = 4 * (floor(Int, T_grad / dwell) ÷ 4)
adc = ADC(Nadc, (Nadc - 1) * dwell, D + dwell / 2);

# The block holds the leading dead time, the whole acquisition window, and the trailing
# dead time, rounded up to the block duration raster.
T_block = ceil((D + Nadc * dwell + lim.ADC_dead_time) / lim.DUR_Δt) * lim.DUR_Δt;

# `Grad(A, T, rise, fall, delay)` with a vector `A` plays a uniformly sampled waveform.
gx = Grad(result.waveform.g[:, 1], T_grad, 0, 0, D)
gy = Grad(result.waveform.g[:, 2], T_grad, 0, 0, D)
gz = Grad(0, 0)
readout = Sequence([gx; gy; gz;;], [RF(0, 0);;], [adc], [T_block])

seq = Sequence(scanner)
@addblock seq += readout

p6 = plot_seq(seq; height=320)
#jl display(p6);

# Koma computes the k-space trajectory from the gradients themselves, so this is an
# independent check that the conversion kept the spiral intact.
p7 = KomaMRI.plot_kspace(seq; view_2d=true, width=400, height=400)
#jl display(p7);

# ## Exporting to Pulseq
#
# Finally, `write_seq` checks the timing and hardware limits against `scanner` and writes
# the readout to a Pulseq file.
seq.DEF["Name"] = "optiks_spiral"
seq_file = joinpath(tempdir(), "optiks_spiral.seq")
write_seq(seq, seq_file; sys=scanner)

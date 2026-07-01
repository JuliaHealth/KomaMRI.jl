# # Building and Exporting a Pulseq GRE Sequence
#
# Pulseq examples usually build sequences from `mr.make*` event constructors.
# KomaMRI now has matching PulseDesigner constructors, but the resulting events
# are native Koma objects: they can be plotted, exported to Pulseq, or used in
# Koma simulation workflows.
#
# This tutorial follows Pulseq's compact
# [`mini_gre.m`](https://github.com/pulseq-admin/pulseq/blob/master/matlab/demoUnsorted/mini_gre.m):
# a slice-selective RF pulse, phase prephasing, Cartesian readout, and a small
# phase-encoding loop.

using KomaMRI
using KomaMRI.PulseDesigner # May move to its own package.
using Unitful # Lets us write values with units like 20u"mT/m".

# We start with scanner limits and sequence parameters.
sys = Scanner(
    B1=50u"μT", Gmax=20u"mT/m", Smax=100u"T/m/s", ADC_Δt=100u"ns",
    RF_ring_down_time=20u"μs", RF_dead_time=100u"μs", ADC_dead_time=10u"μs",
);

# `Scanner` keyword arguments are converted independently. Plain numbers are SI;
# Unitful quantities can be mixed in per field when the physical unit is clear.
FOV = 256u"mm"
slice_thickness = 5u"mm"
Nx = Ny = 128
flip_angle = 15u"deg"
TE = 8u"ms"
TR = 22u"ms"
rf_duration = 4u"ms"
readout_time = 6.4u"ms"
prephaser_time = 2u"ms";

# ## Pulseq-style constructors
#
# Pulseq's `seq.addBlock(...)` appends one block: a set of events played
# together. A Koma `Sequence` is the ordered list of those blocks.
#
# - `make_*` creates events or event tuples, like Pulseq's `mr.make*`.
# - `build_*` creates a `Sequence` with the event already placed into block(s).
#
# We use `make_*` throughout so the blocks below follow the same structure as
# the Pulseq mini GRE example.

rf, gz, gz_reph = make_sinc_pulse(
    flip_angle; duration=rf_duration, slice_thickness, apodization=0.5,
    time_bw_product=4, sys, use=Excitation(),
);

Δk = 1 / FOV
readout_area = Nx * Δk

gx = make_trapezoid(; flat_area=readout_area, flat_time=readout_time, sys)
adc = make_adc(Nx; duration=gx.T, delay=gx.rise, sys)
gx_area = area(gx) * u"T*s/m"
gx_pre = make_trapezoid(; area=-gx_area / 2, duration=prephaser_time, sys);

gy_pre = make_trapezoid(; area=Ny / 2 * Δk, duration=prephaser_time, sys)
phase_scales = (0:(Ny - 1)) ./ (Ny / 2) .- 1;

delay_te = make_delay(
    round_to_raster(
        to_SI(TE) - (gz.T / 2 + gz.fall + dur(gx_pre) + dur(gx) / 2),
        sys.GR_Δt,
    ),
)
delay_tr = make_delay(
    round_to_raster(
        to_SI(TR) - (dur(gx_pre) + dur(gz) + dur(gx) + dur(delay_te)),
        sys.GR_Δt,
    ),
);

# ## Building the sequence
#
# `@addblocks` appends each statement as a sequence block. Each loop iteration
# is one TR.
function mini_gre_sequence()
    seq = Sequence(sys)
    @addblocks for pe_scale in phase_scales
        seq += (rf, z=gz)
        seq += (x=gx_pre, y=pe_scale * gy_pre, z=gz_reph)
        seq += delay_te
        seq += (x=gx, adc)
        seq += delay_tr
    end
    return seq
end

seq = mini_gre_sequence()

p1 = plot_seq(seq; range=[0, to_SI(TR) * 1e3], slider=true, height=320)
#jl display(p1);

# The sequence is also a normal Koma sequence, so we can inspect its k-space
# trajectory directly.
p2 = plot_kspace(seq; view_2d=true, width=400, height=400)
#jl display(p2);

# ## Exporting to Pulseq
#
# The same object can be written as a `.seq` file. As in Pulseq, sequence
# definitions can store metadata such as FOV and sequence name.
# `write_seq` checks timing and hardware limits by default. Passing `sys` makes
# those checks and the Pulseq raster use the scanner above instead of `seq.DEF`.
seq.DEF["FOV"] = [FOV, FOV, slice_thickness] .|> to_SI
seq.DEF["Name"] = "mini_gre"
seq_file = joinpath(tempdir(), "mini_gre.seq")
write_seq(seq, seq_file; sys)

# The exported sequence matches MATLAB Pulseq's `mini_gre.m` output at
# floating-point precision.
# To compare exported files against MATLAB Pulseq's default text precision, use:
# `write_seq(seq, seq_file; significant_digits=6, shape_significant_digits=9)`.

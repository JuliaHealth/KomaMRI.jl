# # Building and Exporting a Pulseq GRE Sequence
#
# Pulseq examples usually build sequences from `mr.make*` event constructors.
# KomaMRI now has matching PulseDesigner constructors, but the resulting events
# are native Koma objects: they can be plotted, exported to Pulseq, or used in
# Koma simulation workflows.
#
# This tutorial follows the structure of Pulseq's compact `mini_gre.m` example:
# a slice-selective RF pulse, phase prephasing, Cartesian readout, and a small
# phase-encoding loop.

using KomaMRI
using KomaMRI.PulseDesigner # May move to its own package.
using Unitful # Optional; without units, constructors use SI values.

# We start with scanner limits and sequence parameters.
sys = Scanner(
    Gmax=20u"mT/m", Smax=40u"T/m/s",
    RF_ring_down_time=20u"μs", RF_dead_time=100u"μs", ADC_dead_time=10u"μs",
);

FOV = 256u"mm"
slice_thickness = 5u"mm"
Nx = Ny = 32
flip_angle = 15u"deg"
TE = 9u"ms"
TR = 22u"ms"
rf_duration = 4u"ms"
readout_time = 6.4u"ms"
prephaser_time = 2u"ms";

# ## Pulseq-style constructors
#
# PulseDesigner follows a simple convention:
#
# - `make_*` creates event objects, such as `Grad`, `RF`, or `ADC`.
# - `build_*` creates a one-block `Sequence`, or a short sequence when the
#   Pulseq event naturally expands to multiple blocks.
#
# For the slice-selective RF pulse, `build_sinc_pulse` creates the RF block and
# the slice rephaser block.

excitation = build_sinc_pulse(flip_angle; duration=rf_duration, slice_thickness, sys, use=Excitation());

γ_unit = γ * u"Hz/T"
Δk = 1 / FOV
readout_area = Nx * Δk / γ_unit

gx = make_trapezoid(; flat_area=readout_area, flat_time=readout_time, sys)
adc = make_adc(Nx; duration=readout_time, delay=gx.rise * u"s", sys)
gx_pre = make_trapezoid(; area=-area(gx) / 2 * u"T*s/m", duration=prephaser_time, sys);

phase_areas = ((0:(Ny-1)) .- (Ny - 1) / 2) .* Δk ./ γ_unit
max_grad_area = maximum(abs.(phase_areas))
gy_pre = make_trapezoid(; area=max_grad_area, duration=prephaser_time, sys);
phase_scales = phase_areas ./ max_grad_area;

# `build_delay` is the sequence-returning version of `make_delay`; it is useful
# when the delay is meant to be appended as its own block.
rf_center = excitation.RF[1].delay + excitation.RF[1].center # Sequence blocks are indexable.
delay_te = build_delay(TE - (dur(excitation) - rf_center + dur(gx_pre) + dur(gx) / 2) * u"s"; sys)
delay_tr = build_delay(TR - (dur(excitation) + dur(gx_pre) + dur(delay_te) + dur(gx)) * u"s"; sys);

# ## Building the sequence
#
# `@addblocks` appends each statement as a sequence block. Each loop iteration
# is one TR.
function mini_gre_sequence()
    seq = Sequence(sys)
    @addblocks for pe_scale in phase_scales
        seq += excitation
        seq += (x=gx_pre, y=pe_scale * gy_pre) # Pulseq: mr.scaleGrad(gyPre, peScale)
        seq += delay_te
        seq += (x=gx, adc)
        seq += delay_tr
    end
    return seq
end

seq = mini_gre_sequence()

p1 = plot_seq(seq; range=[0, TR / u"ms"], slider=true, height=320)
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
seq.DEF["FOV"] = Float64.(ustrip.(u"m", [FOV, FOV, slice_thickness]))
seq.DEF["Name"] = "mini_gre"
seq_file = joinpath(tempdir(), "mini_gre.seq")
write_seq(seq, seq_file; sys)

# Reading it back gives a normal Koma sequence again.
seq_roundtrip = read_seq(seq_file; verbose=false)
seq_roundtrip ≈ seq

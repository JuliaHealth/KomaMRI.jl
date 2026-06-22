using KomaMRI
using KomaMRI.PulseDesigner # May move to its own package.
using Unitful # Optional; without units, constructors use SI values.

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

rf_center = excitation.RF[1].delay + excitation.RF[1].center # Sequence blocks are indexable.
delay_te = build_delay(TE - (dur(excitation) - rf_center + dur(gx_pre) + dur(gx) / 2) * u"s"; sys)
delay_tr = build_delay(TR - (dur(excitation) + dur(gx_pre) + dur(delay_te) + dur(gx)) * u"s"; sys);

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
display(p1);

p2 = plot_kspace(seq; view_2d=true, width=400, height=400)
display(p2);

seq.DEF["FOV"] = Float64.(ustrip.(u"m", [FOV, FOV, slice_thickness]))
seq.DEF["Name"] = "mini_gre"
seq_file = joinpath(tempdir(), "mini_gre.seq")
write_seq(seq, seq_file; sys)

seq_roundtrip = read_seq(seq_file; verbose=false)
seq_roundtrip ≈ seq

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

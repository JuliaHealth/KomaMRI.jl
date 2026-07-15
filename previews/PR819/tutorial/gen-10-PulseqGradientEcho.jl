using KomaMRI
using KomaMRI.PulseDesigner # May move to its own package.
using Unitful # Lets us write values with units like 20u"mT/m".

sys = Scanner(
    B1=50u"μT", Gmax=20u"mT/m", Smax=100u"T/m/s", ADC_Δt=100u"ns",
    RF_ring_down_time=20u"μs", RF_dead_time=100u"μs", ADC_dead_time=10u"μs",
);

FOV = 256u"mm"
slice_thickness = 5u"mm"
Nx = Ny = 128
flip_angle = 15u"deg"
TE = 8u"ms"
TR = 22u"ms"
rf_duration = 4u"ms"
readout_time = 6.4u"ms"
prephaser_time = 2u"ms";

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

seq = Sequence(sys)
@addblock for pe_scale in phase_scales
    seq += (rf, z=gz)
    seq += (x=gx_pre, y=pe_scale * gy_pre, z=gz_reph)
    seq += delay_te
    seq += (x=gx, adc)
    seq += delay_tr
end

p1 = plot_seq(seq; range=[0, to_SI(TR) * 1e3], slider=true, height=320)
display(p1);

p2 = plot_kspace(seq; view_2d=true, width=400, height=400)
display(p2);

seq.DEF["FOV"] = [FOV, FOV, slice_thickness] .|> to_SI
seq.DEF["Name"] = "mini_gre"
seq_file = joinpath(tempdir(), "mini_gre.seq")
write_seq(seq, seq_file; sys)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

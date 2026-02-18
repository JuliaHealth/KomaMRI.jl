using KomaMRI, Suppressor
sys = Scanner(); # default hardware definition

ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF,durRF);

nADC = 8192         # number of acquisition samples
durADC = 250e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay);

seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
p1 = plot_seq(seq; slider=false, height=300)
display(p1);

obj = Phantom(x=[0.], T1=[1000e-3], T2=[100e-3]);

raw = @suppress simulate(obj, seq, sys);

p2 = plot_signal(raw; slider=false, height=300)
display(p2);

obj = Phantom(x=[0.], T1=[1000e-3], T2=[100e-3], Δw=[-2π*100])# and simulate again.

raw = @suppress simulate(obj, seq, sys)
p3 = plot_signal(raw; slider=false, height=300)
display(p3);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

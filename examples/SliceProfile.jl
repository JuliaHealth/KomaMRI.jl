using Koma, Plots, LaTeXStrings
pgfplotsx()
rf_wf = [
    #Waveform am_sg_100_100_0
    0, 413, 859, 1337, 1847, 2389, 2962, 3568, 4204, 4870,
    5565, 6288, 7038, 7814, 8613, 9435, 10277, 11138, 12015, 12906,
    13809, 14721, 15640, 16563, 17487, 18409, 19327, 20237, 21136, 22023,
    22892, 23742, 24570, 25372, 26146, 26889, 27598, 28271, 28905, 29497,
    30046, 30549, 31005, 31411, 31767, 32070, 32320, 32515, 32655, 32739,
    32767, 32739, 32655, 32515, 32320, 32070, 31767, 31411, 31005, 30549,
    30046, 29497, 28905, 28271, 27598, 26889, 26146, 25372, 24570, 23742,
    22892, 22023, 21136, 20237, 19327, 18409, 17487, 16563, 15640, 14721,
    13809, 12906, 12015, 11138, 10277, 9435, 8613, 7814, 7038, 6288,
    5565, 4870, 4204, 3568, 2962, 2389, 1847, 1337, 859, 413,
    0
] ./ (2^16/2 - 1) #signed int 16 bits
Trf = 0.486e-3
Arf = 21e-6 * rf_wf
Gx1 = 12.89e-3
ζ11 = 0.29e-3
ζ21 = 0.37e-3
Tg1 = Trf
delay = Tg1 + ζ11 - Trf
ex1 = Sequence([Grad(0,0); Grad(0,0); Grad(Gx1,Tg1,ζ11,ζ21,0);;],[RF(Arf,Trf,-5000,delay);;],[ADC(0,0)])
ex2 = Sequence([Grad(0,0); Grad(0,0); Grad(Gx1,Tg1,ζ11,ζ21,0);;],[RF(Arf,Trf,    0,delay);;],[ADC(0,0)])
ex3 = Sequence([Grad(0,0); Grad(0,0); Grad(Gx1,Tg1,ζ11,ζ21,0);;],[RF(Arf,Trf, 5000,delay);;],[ADC(0,0)])
Gx2 = -.83e-3
ζ12 = 0.026e-3
ζ22 = 0.28e-3
Tg2 = 5.75e-3
ref = Sequence([Grad(0,0); Grad(0,0); Grad(Gx2,Tg2,ζ12,ζ22,0);;])
#Different slices
seq1 = ex1 + ref
seq2 = ex2 + ref
seq3 = ex3 + ref
#Focus Time
t_focus_aux = (abs(Gx2)*(2Tg2+ζ12+ζ22)-ζ21*Gx1)/(2Gx1)
t_focus = (Tg1+ζ11-t_focus_aux)
#SIM
z = range(-2e-2,2e-2,200)
M1 = simulate_slice_profile(seq1;z)
M2 = simulate_slice_profile(seq2;z)
M3 = simulate_slice_profile(seq3;z)
#Plots
plot_seq(seq2; width=480, height=350, slider=false)
f = γ*Gx1*z*1e-3
plot(f,  abs.(M1.xy),label=L"|M_{xy}| @ \Delta f = -5\, \mathrm{kHz}",
        xlabel=L"\mathrm{kHz}",line=3,title="α=$(KomaMRI.get_flip_angle(seq2.RF[1]))")
plot!(f, abs.(M2.xy),label=L"|M_{xy}| @ \Delta f = +0\, \mathrm{kHz}",line=3)
plot!(f, abs.(M3.xy),label=L"|M_{xy}| @ \Delta f = +5\, \mathrm{kHz}",line=3)

#GIF
# anim = @animate for tf = range(0,7.5e-3,60)
#     simParams = Dict("return_Mag"=>true, "end_sim_at"=>tf, "gpu"=>false)
#     M = simulate(phantom, seq2, sys; simParams)
#     #Plots
#     plot(xx*1e2, abs.(M.xy),label="Mxy",line=4,ylim=(-1,1),title="Simulation ended at $(round(simParams["end_sim_at"]*1e3,digits=3)) ms")
#     plot!(xx*1e2, abs.(M.z),label="Mz",line=:dash)
#     plot!(xx*1e2, real.(M.xy),label="Mx",line=:dash)
#     plot!(xx*1e2, imag.(M.xy),label="My",line=:dash)
# end
# gif(anim, "AnimSliceExc.gif", fps = 10)

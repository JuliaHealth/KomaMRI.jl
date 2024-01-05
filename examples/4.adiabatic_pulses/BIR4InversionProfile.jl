# This document replicates the results of Figure 15g of the paper "The Return of the Frequency Sweep:
# Designing Adiabatic Pulses for Contemporary NMR" by Michael Garwood and Lance DelaBarre.
using KomaMRI, MAT, PlotlyJS, LinearAlgebra, ProgressMeter

RF_wf = matread("./examples/4.adiabatic_pulses/BIR4.mat")
R = 200 #Product duration and bandwidth
B1 = 120e-6; w1_2pi_Hz = γ*B1
Trfs = [2] * 1e-3 #(4:0.1:7.5) * 1e-3 #ms
f0s = R ./ (2Trfs) # R = 2 fmax * Trf
ΔB1s = 0.8:0.2:1.2 # %
#Init
score = zeros(length(f0s), length(Trfs), length(ΔB1s))
p = [scatter() for i=1:length(f0s), j=1:length(Trfs), k=1:length(ΔB1s), l=1:2]
N = prod(size(score))
count = 1
for (i, f0) = enumerate(f0s), (j, Trf) = enumerate(Trfs), (k, ΔB1) = enumerate(ΔB1s)
    b1 = 120e-6 * ΔB1
    B1 = b1 * RF_wf["b1"][:]
    Δf = RF_wf["df"][:] * f0
    dt = Trf / length(B1)
    fmax_sim = 10e3
    Gz = fmax_sim / γ
    seq = Sequence(
        [Grad(0.,0.); Grad(0.,0.); Grad(Gz,Trf,0);;],
        [RF(B1,Trf,Δf,0);;]
        )
    sim_params = Dict{String,Any}("Δt_rf"=>dt)

    z = range(-1, 1, 200)
    M = simulate_slice_profile(seq; sim_params, z)

    f = γ*Gz*z
    p[i,j,k,1] = scatter(x=f,y=(1 .- M.z)/2,name="B1/B1_nom=$ΔB1")
    # p[i,j,k,2] = scatter(x=f,y=abs.(M.xy),name="Mxy, B1/B1_nom=$ΔB1")
    println("################ $count / $N = $(round(count/N*100; digits=3)) % ################")
    global count = count + 1
end
# Heatmap
plot(p[:], Layout(title="Fraction Refocused (BIR-4, Trf=$(Trfs[1]*1e3)ms, TBP=Trf*BW=$(R))"))

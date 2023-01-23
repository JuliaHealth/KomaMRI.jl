# This document replicates the results of Figure 15g of the paper "The Return of the Frequency Sweep: 
# Designing Adiabatic Pulses for Contemporary NMR" by Michael Garwood and Lance DelaBarre.
using KomaMRI, MAT, PlotlyJS, LinearAlgebra, ProgressMeter

RF_wf = matread("./examples/4.adiabatic_pulses/BIR4.mat")
for R = [50] #Product duration and bandwidth
    B1_nom = 13.5e-6; w1_2pi_Hz = γ*B1_nom
    Trfs = [7] * 1e-3  # (4:0.1:7.5) * 1e-3 #ms
    f0s = R ./ (2Trfs) # R = (2 fmax) * Trf
    ΔB1s = Array(0.3:0.01:1.5) # %
    NB0s = 200
    B0_max = 100 #+-100Hz
    Gz = B0_max / γ
    z = range(-1, 1, NB0s)
    B0s = Array(γ*Gz*z)

    #Init
    score = zeros(length(f0s), length(Trfs), length(ΔB1s))
    MagXY = zeros(ComplexF64, NB0s, length(ΔB1s))
    MagZ = zeros(ComplexF64, NB0s, length(ΔB1s))
    N = prod(size(score))
    counter = 1
    for (i, f0) = enumerate(f0s), (j, Trf) = enumerate(Trfs), (k, ΔB1) = enumerate(ΔB1s)
        b1 = B1_nom * ΔB1
        B1 = b1 * RF_wf["b1"][:]
        Δf = RF_wf["df"][:] * f0
        dt = Trf / length(B1)
        T90 = .5e-3
        B190 = 90 / (360 * γ * T90)
        rf90 = Sequence(
            [Grad(0.,0.); Grad(0.,0.); Grad(Gz,T90,0);;], 
            [RF(ΔB1*B190,T90,0,0);;]
            )
        bir4 = Sequence(
            [Grad(0.,0.); Grad(0.,0.); Grad(Gz,Trf,0);;], 
            [RF(B1,Trf,Δf,0);;]
            )
        seq = Sequence()
        # seq += rf90
        seq += bir4
        # seq += rf90
        simParams = Dict{String,Any}("Δt_rf"=>dt)

        M = simulate_slice_profile(seq; simParams, z)
        # display(plot_seq(seq))

        MagXY[:, counter] = M.xy
        MagZ[:, counter] = M.z
        println("################ $counter / $N = $(round(counter/N*100; digits=3)) % ################")
        counter = counter + 1
    end
    # Heatmap
    fmax = round(f0s[1] * 1e-3; digits=2)
    p1 = plot(heatmap(y=B0s, x=ΔB1s, z=real.(MagXY), colorscale="Jet"), Layout(title="MX BIR-4 Trf=7 ms R=$R fmax=$fmax kHz"))
    p2 = plot(heatmap(y=B0s, x=ΔB1s, z=imag.(MagXY), colorscale="Jet"), Layout(title="MY BIR-4 Trf=7 ms R=$R fmax=$fmax kHz"))
    p3 = plot(heatmap(y=B0s, x=ΔB1s, z=real.(MagZ), colorscale="Jet"), Layout(title="MZ BIR-4 Trf=7 ms R=$R fmax=$fmax kHz"))
    display([p1; p2; p3])
end
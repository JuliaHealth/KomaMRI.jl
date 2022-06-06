# Experiment to see how the motion-compensated diffusion sequences
# affects microstructural signals.

using Koma, JLD2, ProgressMeter
# using PlotlyJS: savefig

Ls = (1:4:120)*1e-6
bs = 0:.25:10
E = zeros(ComplexF64,length(Ls),length(bs))
Emc = zeros(ComplexF64,length(Ls),length(bs))
@load "M0MC.jld2"
DIF1 = DIF
@load "M2MXMC.jld2"
DIF2 = DIF
##Timings
@showprogress for (j,b) in enumerate(bs)
    for (i,L) in enumerate(Ls)
        μ = KomaMRI.Planes(L)
        b1 = KomaMRI.get_bvalue(DIF1)*1e-6
        seq = DIF1 * b
        E[i,j] = KomaMRI.SignalE(μ, seq)
        b2 = KomaMRI.get_bvalue(DIF2)*1e-6
        seq = DIF2 * sqrt(b1/b2) * b
        Emc[i,j] = KomaMRI.SignalE(μ, seq)
    end
end
@save "ResultsMC.jld2" E Emc
##
# p = KomaMRI.plot_grads_moments(DIF1,title="M0, b=$(round(KomaMRI.get_bvalue(DIF1)*1e-6, digits=2)) s/mm2",mode="hd")
# PlotlyJS.savefig(p,"M0.pdf";width=900,height=250)
# ##
# b1 = KomaMRI.get_bvalue(DIF1)*1e-6
# b2 = KomaMRI.get_bvalue(DIF2)*1e-6
# DIF = DIF2 * sqrt(b1/b2)
# p = KomaMRI.plot_grads_moments(DIF,title="M2 M1 MX MC, b=$(round(KomaMRI.get_bvalue(DIF)*1e-6, digits=2)) s/mm2",mode="hd")
# PlotlyJS.savefig(p,"M1M2MC.pdf";width=900,height=250)
##
@load "ResultsMC.jld2" E Emc

using Plots#, Interact, Blink
b1 = KomaMRI.get_bvalue(DIF1)*1e-6
# w = Window()
anim = @animate for l = 1:9
    D = 2e-9
    b = b1 * 1e6 * bs[l]
    p = plot(Ls*1e6,abs.(E[:,l]),label="M0",legend=:topright,ylim=(0,1),xlim=(0,30))
    plot!(Ls*1e6,abs.(Emc[:,l]),label="M2-Maxwell")
    plot!(Ls*1e6,abs.(ones(size(Emc[:,l])) .* exp(-b*D)),label="Free diffusion")
    xlabel!("L [μm]")
    ylabel!("E(q)")
    title!("Comparison (b=$(round(b*1e-6)) s/mm2)")
end
# body!(w,p)
gif(anim, "MC_L.gif", fps = 1)
# Plots.savefig(p,"Eq.pdf")

## Plot
anim = @animate for l = 1:10
    D = 2e-9
    p = plot(b1*bs, abs.(E[l,:]),label="M0",legend=:topright,xlim=(0,2000))# ,yaxis=:log)
    plot!(b1*bs, abs.(Emc[l,:]),label="M2-Maxwell")
    plot!(b1*bs, exp.(-b1*bs*1e6*D),label="Free diffusion")
    xlabel!("b [s/mm2]")
    ylabel!("E(q)")
    title!("Comparison (L=$(round(Ls[l]*1e6)) um)")
end
gif(anim, "MC_B.gif", fps = 2)

## Effect of rotationa angle
θs = -π/2:.1:π/2
E = zeros(ComplexF64,length(θs))
Emc = zeros(ComplexF64,length(θs))
f = 2
L = 20e-6
@showprogress for (i,θ) in enumerate(θs)
        μ = KomaMRI.Planes(L)
        R = KomaMRI.rotz(θ)
        b1 = KomaMRI.get_bvalue(DIF1)*1e-6
        seq = DIF1 * f * R
        E[i] = KomaMRI.SignalE(μ, seq)
        b2 = KomaMRI.get_bvalue(DIF2)*1e-6
        seq = DIF2 * sqrt(b1/b2) * f * R
        Emc[i] = KomaMRI.SignalE(μ, seq)
end
##
using Plots
f = 2
D = 2e-9
b = KomaMRI.get_bvalue(DIF1)
p = plot(θs,abs.(E[:]),label="M0",legend=:topright,ylim=(0,1))
plot!(θs,abs.(Emc[:]),label="M2-Maxwell")
plot!(θs,abs.(ones(size(Emc[:])) .* exp(-b*f*D)),label="Free diffusion")
xlabel!("θ [rad]")
ylabel!("E(q)")
title!("Comparison (L=20 um, b=$(round(b*1e-6*f)) s/mm2)")
using MRIsim, Plots
using MRIsim: γ, Planes, SignalE, DIF_base, plot_grads

#Constants
𝒊 = 1im; D = 2e-9 # m2/s
#DIF Sequence
T = 80e-3; δ = 10e-3;
DIF = DIF_base(1,T-δ,δ); #plot_grads(DIF)
#Restricted planes with Short Pulse Approximation (SPA) - from my thesis 
Erest(L,q) = sum([(n==0 ? 1 : 2)*exp(-D*(T-δ)*(n*π/L)^2)*sinc.(L*q .+ n/2).^2 ./(1 .- n./(2*L*q)).^2 for n = 0:100])
#Free diffusion - SPA
Efree(q) = exp.(-4π^2*D*q^2)

## Comparison
gr()
plot()
g = (0:1000)*1e-3
q = γ*g*δ #q-values
color = [:blue, :red, :green]
for (i, L) = enumerate([2, 7, 90]*1e-6)
    Ll = round(Int,L*1e6)
    μ = Planes(L,D)
    #Laplacian Eigen Functions for different gradient strengths
    E = [SignalE(μ, gi*DIF) for gi ∈ g] 
    plot!(q*1e-6, abs.(E), linecolor=color[i], label="L = "*string(Ll)*" μm", yaxis=:log, legend=:outertopright)
    plot!(q*1e-6, abs.(Erest.(L,q)), linecolor=color[i], label="L = "*string(Ll)*" μm", linestyle=:dash,yaxis=:log, legend=:outertopright)
end
plot!(q*1e-6, abs.(Efree.(q)),linecolor=:purple,label="L = ∞ μm",yaxis=:log, legend=:outertopright)
ylims!(1e-6,1)
ylabel!("E(q)")
xlabel!("q [μm⁻¹]")
title!("Laplacian Eigen Functions vs Short Pulse Approximation")
#savefig("./src/examples/Figures/PGSEvsSPA.pdf")
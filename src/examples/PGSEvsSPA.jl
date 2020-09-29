using MRIsim, Plots
using MRIsim: γ, Planes

Λ, B, _, U = Planes(30);
𝒊 = 1im; D = 2e-9 # m2/s
T = 80e-3; δ = 10e-3;
#Pulse Gradient Spin Echo (PGSE)
Eq(p,q) = (exp(-(p*Λ .+ 𝒊*q*B)*δ/T)*exp(-p*Λ*(T-2δ)/T)*exp(-(p*Λ .- 𝒊*q*B)*δ/T))[1,1]
#Short Pulse Approximation (SPA)
Et(L,q) = sum([(n==0 ? 1 : 2)*exp(-D*(T-δ)*(n*π/L)^2)*sinc.(L*q .+ n/2).^2 ./(1 .- n./(2*L*q)).^2 for n = 0:100])
#Free diffusion - SPA
Efree(q) = exp.(-4*π^2*D*q^2)
## Comparison
pgfplotsx()
plot()
g = (0:1000)*1e-3
q = γ*g*δ
color = [:blue, :red, :green]
for (i, L) = enumerate([2, 7, 90]*1e-6)
    p = D*T/L^2
    println(round(p,digits=2))
    qun = 2π*γ*g*T*L
    Ll = round(Int,L*1e6)
    plot!(q*1e-6, abs.(Eq.(p,qun)), linecolor=color[i], label=L"L=\,"*string(Ll)*L"\,\mu"*"m", yaxis=:log, legend=:outertopright)
    plot!(q*1e-6, abs.(Et.(L,q)), linecolor=color[i], label=L"L=\,"*string(Ll)*L"\,\mu"*"m", linestyle=:dash,yaxis=:log, legend=:outertopright)
end
plot!(q*1e-6, abs.(Efree.(q)),linecolor=:purple,label=L"L=\infty\,\mu"*"m",yaxis=:log, legend=:outertopright)
ylims!(1e-6,1)
ylabel!(L"E(q)")
xlabel!(L"q\,[\mu \mathrm{m}^{-1}]")
title!("Laplacian Eigen Functions vs Short Pulse Approximation")
#savefig("./src/examples/Figures/PGSEvsSPA.pdf")
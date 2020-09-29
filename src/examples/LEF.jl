using MRIsim, Plots
using MRIsim: γ, Planes, DIF_base, plot_grads, Grad_fun, dur, gpu
## PGSE
𝒊 = 1im; D = 2e-9 # m2/s
T = 80e-3; δ = 30e-3;
Gmax = 30e-3;
DIF = DIF_base(Gmax,T-δ,δ)
plot_grads(DIF)
## Signal decay
Λ, B = Planes(30);
L = 20e-6
p = D*T/L^2
q = 2π*γ*T*L
Ex = *([exp(-(p*Λ .+ 𝒊*q*g.A*B)*g.T/T) for g = DIF.GR[1,:]]...)[1,1]
Ey = *([exp(-(p*Λ .+ 𝒊*q*g.A*B)*g.T/T) for g = DIF.GR[2,:]]...)[1,1]
E = Ex*Ey
## FWF
fwf = Grad_fun(x-> Gmax*sin(2π*x/T), T, 600)
DIF = Sequence([copy(fwf); 0*copy(fwf)])
τ = dur(DIF)
plot_grads(DIF)
## Signal decay
Λ, B = Planes(30);
Ex = *([exp(-(p*Λ .+ 1im*q*g.A*B)*g.T/τ) for g = DIF.GR[1,:]]...)[1,1]
Ey = *([exp(-(p*Λ .+ 𝒊*q*g.A*B)*g.T/τ) for g = DIF.GR[2,:]]...)[1,1]
E = Ex*Ey


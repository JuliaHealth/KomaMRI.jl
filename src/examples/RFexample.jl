using MRIsim
using MRIsim: γ, Q, RF_fun, get_grads
# Example
G = 30e-3
T = 3e-3; N = 200
B1 = 15e-6;

Rfs = RF_fun(t->B1*sinc.(8(t-T/2)/T),T)
SEQ = Sequence([Grad(G,T);Grad(0,T)],Rfs)

# TODO
function get_RFs(obj::Phantom,seq::Sequence,t)
    Δt = t[2]-t[1]
    Gx, Gy = get_grads(seq,t)
    φ = -2π*γ*Δt*sqrt.(abs.(B1e).^2) ; φ[φ.==0] .= 1e-17
    n =  2π*γ*Δt*[[real(B1e[i]) imag(B1e[i]) obj.x]./abs(φ[i]) for i = 1:N]
    [Q(φ[i],n[i]) for i=1:N]
end

Qs = [Q(φ[i],n[i]) for i=1:N]
Qt = *(Qs...) #Total rotation matrix

Mxy, Mz = 0. + 0. *im, 1. #[0,0,1]
M = [Mxy, Mz]
Mt = Qt*M
Mz = [(*(Qs[1:i]...)*M)[2] for i=2:N]
Mxy = [(*(Qs[1:i]...)*M)[1] for i=2:N]

#Sphere
u = range(0,stop=2π,length=50);
v = range(0,stop=π,length=50);
x = .99 * cos.(u) * sin.(v)';
y = .99 * sin.(u) * sin.(v)';
z = .99 * ones(100) * cos.(v)';

## Plots
using Plots, LaTeXStrings
pgfplotsx()
plot(real.(Mxy),imag.(Mxy),real.(Mz),linewidth=5,
label=L"M(t)",xlabel="x",ylabel="y",zlabel="z")
scatter!([real.(Mxy[1])],[imag.(Mxy[1])],[real.(Mz[1])],
          marker=:xcross,markersize=1.5,markerstrokecolor=:auto,label="Start")
scatter!([real.(Mxy[end])],[imag.(Mxy[end])],[real.(Mz[end])],
          marker=:xcross,markersize=1.5,markerstrokecolor=:auto,label="End")
wireframe!(x,y,z,color=:gray)
##
using Plots
pgfplotsx()
plot()

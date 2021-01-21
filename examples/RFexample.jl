using MRIsim
using MRIsim: γ, Q, RF_fun, get_grads
# Example
G = 30e-3 # 30 mT/m
T = 3e-3; # 3 ms
B1 = 15e-6; # 15 μT

RFs = RF_fun(t->B1*sinc.(8(t-T/2)/T),T) #sinc pulse of duration T [s]
B1e = getproperty.(SEQ.RF,:A)
ΔT = getproperty.(SEQ.RF,:T)
Bx, By = real.(B1e), imag.(B1e)
Bz = zeros(size(B1e)) # Bz = 0 in the rotating frame implies a non-selective excitation
# In genral, Bz = Δω(x)/γ + G ⋅ x

φ = -2π*γ * ΔT .* sqrt.(abs.(B1e).^2 .+ abs.(Bz).^2) # angle of rotation 
φ[φ.==0] .= 1e-17; # removes problems when dividing by φ
n =  2π*γ * ΔT .* [Bx By Bz]./abs.(φ) # axis of rotation

Qs = [Q(φ[i],n[i,:]) for i=1:N] # hard-pulse approximation for every RF element in RFs
Qt = *(Qs...) # Total rotation matrix

Mxy, Mz = 0. + 0. *im, 1. # [0,0,1]
M0 = [Mxy, Mz] # Initial magnetization
Mt = Qt*M0 # Final magnetization
Mz = [(*(Qs[1:i]...)*M0)[2] for i=2:N]
Mxy = [(*(Qs[1:i]...)*M0)[1] for i=2:N]

## Sphere
u = range(0,stop=2π,length=50);
v = range(0,stop=π,length=50);
x = .99 * cos.(u) * sin.(v)';
y = .99 * sin.(u) * sin.(v)';
z = .99 * ones(100) * cos.(v)';

## Plots
using Plots
plotlyjs() #change backend
surface(x[:],y[:],z[:])
plot!(real.(Mxy),imag.(Mxy),real.(Mz),
    linewidth=5,label=L"M(t)",xlabel="x",ylabel="y",zlabel="z",
    xlims=(-1,1),ylims=(-1,1),zlims=(-1,1))
scatter!([real.(Mxy[1])],[imag.(Mxy[1])],[real.(Mz[1])],
          marker=:xcross,markersize=1.5,markerstrokecolor=:auto,label="Start")
scatter!([real.(Mxy[end])],[imag.(Mxy[end])],[real.(Mz[end])],
          marker=:xcross,markersize=1.5,markerstrokecolor=:auto,label="End")


## TODO
# - A function that translates from Phantom and Sequence to a rotation matrix for each spin, something like:
# function get_rf_rotations(obj::Phantom, seq::Sequence, t)
#     Δt = t[2]-t[1]
#     Gx, Gy = get_grads(seq,t)
#     B1e = get_rfs(seq,t) #TODO: uniformly sample RF pulses with respect to the simulation sampling period Δt
#     for spin ∈ obj # this an example, it will not work
#         φ = -2π*γ * Δt .* sqrt.(abs.(B1e).^2 .+ abs.(Bz).^2) # angle of rotation 
#         φ[φ.==0] .= 1e-17; # removes problems when dividing by φ
#         n =  2π*γ * Δt .* [Bx By Bz]./abs.(φ) # axis of rotation
#         Qtotal = prod([Q(φ[i],n[i,:]) for i=1:N]) # hard-pulse approximation for every RF element in RFs
#     end
#     Qtotal
# end
# - Obtaining rotation matrices efficiently: as we would NEVER acquire during an RF pulse we do not need
#   to simulate intermidiate steps. This is due to the fact that we need to change between Rx and Tx with a switch (hardware limitation).
# - I suggest to simulate RF blocks separatly, and then to use run_sim2D_spin() in between.
#   The only thing we need to add is a way to follow the T1 decay during the non-RF parts of the simulation.
#   The magnetization [Mx + iMx, Mz] should decay as follows (the first part is already in run_sim2D_spin()):
#                    [exp(-Δt/T2)(Mx + iMx), exp(-Δt/T1) Mz + ρ (1 - exp(-Δt/T1))].
# - Movement of the phantom during excitation? (using obj.x + obj.ux(obj.x,obj.y))
# - Effects of T1 and T2 during an RF block: For each Δt we do RF and then T1 & T2 decay.
# - Plotting function MRIsim.plot_grads(seq) in Display.jl is not currently plotting RF pulses
# - Pre-defined RF waveforms.
# - Shinnar-Le Roux pulse desing? 

# Precursor of the RF excitation code with feature directions, most of the comments below are already implemented

using KomaMRI
using KomaMRI: γ, Q, RF_fun, get_grads, Un, cross
using Plots, LaTeXStrings
gr(size = (800,800))

# Example
G = 30e-3 # 30 mT/m
T = 3e-3; # 3 ms
B1 = 15e-6; # 15 μT

N = 300
RFs = RF_fun(t->B1*sinc.(4(t-T/2)/T),T,N) #sinc pulse of duration T [s]
B1e = RFs.A
α = π / 2
ΔT = RFs.T
c = α / (2π*γ*sum(B1e.*ΔT))
c = c.re
RFs.A .= RFs.A * c

plot(real.(RFs.A), linewidth=5, dpi=200)

## In general, Bz = Δω(x)/γ + G ⋅ x
B1e = RFs.A
ΔT = RFs.T
Bx, By = real.(B1e), imag.(B1e)

# Plot
df = .2 #.02
ff = -1:df:1
colors = range(HSV(0,1,1), stop=HSV(-360,1,1), length=length(ff)) # inverse rotation

anim = @animate for j ∈ reverse(2:N)
plot()
for (i,f) = enumerate(ff)
Bz = zeros(size(B1e)) .+ f*.01e-3 # Bz = 0 in the rotating frame implies a non-selective excitation

φ = -2π*γ * ΔT .* sqrt.(abs.(B1e).^2 .+ abs.(Bz).^2) # angle of rotation 
φ[φ.==0] .= 1e-17; # removes problems when dividing by φ
n =  2π*γ * ΔT .* [Bx By Bz]./abs.(φ) # axis of rotation
Qs = [Q(φ[i],n[i,:]) for i=1:N] # hard-pulse approximation for every RF element in RFs
Qt = *(Qs...) # Total rotation matrix

Mxy, Mz = 0. + 0. *im, 1. # [0,0,1]
M0 = Mag(Mxy, Mz) # Initial magnetization
Mt = Qt*M0 # Final magnetization
Mz = [(*(Qs[1:i]...)*M0)[2] for i=2:N-j]
Mxy = [(*(Qs[1:i]...)*M0)[1] for i=2:N-j]

# Plots
plot!(real.(Mxy),imag.(Mxy),real.(Mz),
    linewidth=5,xlabel="x",ylabel="y",zlabel="z",
    xlims=(-1,1),ylims=(-1,1),zlims=(-1,1),dpi=200,label=:none,
    linecolor=colors[i],subplot=1)
end

plot!(real.(RFs.A)[1:N-j], linewidth=5, dpi=200,
    xlims=(0,N),ylims=(-.5c*B1,1.1c*B1), ticks = nothing,
    inset_subplots=bbox(0.05, 0.05, 0.5, 0.25, :bottom, :right),
    subplot=2,label=:none)

end
#gif(anim, "anim_fps60.gif", fps = 60)

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
#   to store intermidiate steps. This is due to the fact that we need to change between Rx and Tx with a switch (hardware limitation).
# - I suggest to simulate RF blocks separatly, and then to use run_sim2D_spin() in between.
#   The only thing we need to add is a way to follow the T1 decay during the non-RF parts of the simulation.
#   The magnetization [Mx + iMx, Mz] should decay as follows (the first part is already in run_sim2D_spin()):
#                    [exp(-Δt/T2)(Mx + iMx), exp(-Δt/T1) Mz + ρ (1 - exp(-Δt/T1))].
# - Movement of the phantom during excitation? (using obj.x + obj.ux(obj.x,obj.y))
# - Effects of T1 and T2 during an RF block: For each Δt we do RF and then T1 & T2 decay.
# - Plotting function KomaMRI.plot_grads(seq) in Display.jl is not currently plotting RF pulses
# - Pre-defined RF waveforms.
# - Shinnar-Le Roux pulse desing? 

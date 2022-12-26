#This example replicates the results of page 41 on the book "Handbook of MRI Pulse Sequences" by Bernstein et al.
using KomaMRI, PlotlyJS

sys = Scanner()
B1 = 13e-6
Trf = 3.2e-3
zmax = 2e-2
z = range(-zmax, zmax, 400)
fmax = 15e3
Gz = fmax / (γ * zmax)
seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
α = KomaMRI.get_flip_angles(seq)[1]
α_desired = 120 + 0im
seq = (α_desired / α) * seq #Scaling pusle to have a flip angle of 120
plot_seq(seq; max_rf_samples=Inf)

simParams = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, simParams)

f = γ * Gz * z
s1 = scatter(x=f, y=abs.(M.xy), name="|Mxy1|")
p1 = plot([s1], Layout(title="120 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]"))
#This example replicates the results of page 40 on the book "Handbook of MRI Pulse Sequences" by Bernstein et al.

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
α = KomaMRI.get_flip_angles(seq)[1]
α_desired = 30 + 0im
seq = (α_desired / α) * seq #Scaling pusle to have a flip angle of 120
p3 = plot_seq(seq; max_rf_samples=Inf, slider=false)

simParams = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, simParams)

f = γ * Gz * z
s1 = scatter(x=f, y=real.(M.xy), name="Mx")
s2 = scatter(x=f, y=imag.(M.xy), name="My")
p2 = plot([s1,s2], Layout(title="$(abs(α_desired)) deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]"))
[p3 p2]
# # Sinc Slice Profile

filename = last(splitpath(@__FILE__)) # hide
isFileMD = occursin(".md", filename) # hide
isFileJL = occursin(".jl", filename) # hide

#md # This example replicates the results of page 41 on the book
#md # "Handbook of MRI Pulse Sequences" by Bernstein et al.

#md # First of all, let's use the KomaMRI package and define the default scanner.

using KomaMRI
sys = Scanner() # default hardware definition
nothing # hide

#md # Let's define some general parameters.

B1 = 13e-6
Trf = 3.2e-3
zmax = 2e-2
z = range(-zmax, zmax, 400)
fmax = 15e3
Gz = fmax / (γ * zmax)
f = γ * Gz * z
nothing # hide

#md # Define the sequence with a flip angle of 120°.

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
α = KomaMRI.get_flip_angles(seq)[1]
α_desired = 120 * 1im
seq = (α_desired / α) * seq    # Scaling pusle to have a flip angle of 120
p1 = plot_seq(seq; max_rf_samples=Inf)
if isFileMD savefig(p1, "../assets/41-seq.html") end # hide
if isFileJL display(p1) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/41-seq.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Perform the simulation of a slice.

simParams = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, simParams)

using PlotlyJS # hide
sa = scatter(x=f, y=abs.(M.xy), name="|Mxy1|") # hide
pa = plot([sa], Layout(title="120 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]")) # hide
if isFileMD savefig(pa, "../assets/4a-profile.html") end # hide
if isFileJL display(pa) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/4a-profile.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Define the sequence with a flip angle of 120°.

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
α = KomaMRI.get_flip_angles(seq)[1]
α_desired = 30 + 0im
seq = (α_desired / α) * seq  # Scaling pusle to have a flip angle of 120
p2 = plot_seq(seq; max_rf_samples=Inf, slider=false)
if isFileMD savefig(p2, "../assets/42-seq.html") end # hide
if isFileJL display(p2) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/42-seq.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Perform the simulation of a slice.

simParams = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, simParams)

using PlotlyJS # hide
s1 = scatter(x=f, y=real.(M.xy), name="Mx") # hide
s2 = scatter(x=f, y=imag.(M.xy), name="My") # hide
pb = plot([s1,s2], Layout(title="30 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]")) # hide
if isFileMD savefig(pb, "../assets/4b-profile.html") end # hide
if isFileJL display(pb) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/4b-profile.html" style="width:100%; height:320px;"></object>
#md # ```

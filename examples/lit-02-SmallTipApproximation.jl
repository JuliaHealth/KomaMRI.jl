# # Small Tip Angle Approximation

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](@__REPO_ROOT_URL__/examples/lit-02-SmallTipApproximation.jl)

filename = last(splitpath(@__FILE__)) # hide
isFileMD = occursin(".md", filename) # hide
isFileJL = occursin(".jl", filename) # hide


#md # !!! note
#md #     This example replicates the results of page 41 on the book
#md #     "Handbook of MRI Pulse Sequences" by Bernstein et al.

using KomaMRI # hide
sys = Scanner() # hide
sys.Smax = 50 # hide

#md # In this example we will showcase the small tip angle approximation. 
#md # For this, we will simulate a slice profile for spins' postitions ``z\in[-2,\,2]\,\mathrm{cm}``
#md # and will choose a gradient ``G_{z}`` so their frequencies are mapped to ``f\in[-5,\,5]\,\mathrm{kHZ}``.
#md # To start, we define an RF pulse with a flip angle of 30 deg and pulse duration of ``T_{\mathrm{rf}}=3.2\,\mathrm{ms}``.

B1 = 4.92e-6
Trf = 3.2e-3
zmax = 2e-2
fmax = 5e3
z = range(-zmax, zmax, 400)
Gz = fmax / (γ * zmax)
f = γ * Gz * z
α = KomaMRI.get_flip_angles(seq)[1] #Approx 30 deg

#md # The designed RF pulse will look like this (see figure below), 
#md # where the gradient after the RF pulse has half of the area of the first one to refocuse the spins' phase after the excitation.  

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
p2 = plot_seq(seq; max_rf_samples=Inf, slider=false)
if isFileMD savefig(p2, "../assets/42-seq.html") end # hide
if isFileJL display(p2) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/42-seq.html" style="width:100%; height:380px;"></object>
#md # ```

#md # Now we will perform the simulation using the function [`simulate_slice_profile`](@ref). 
#md # Note that we modified `Δt_rf` in `simParams` to match the resolution of the waveform.

simParams = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, simParams)

# Simulation results
using PlotlyJS # hide
s1 = scatter(x=f, y=real.(M.xy), name="Mx") # hide
s2 = scatter(x=f, y=imag.(M.xy), name="My") # hide
# Small tip approximation
dat = seq.RF.A[1] # hide
N = length(dat) # hide
dat_pad = [zeros(floor(Int64,N)); dat; zeros(floor(Int64,N))] # hide
N_pad = length(dat_pad) # hide
U = 1 / (Trf) * N / N_pad #hide
u = range(0, (N_pad - 1) * U; step=U) # hide
u = u .- maximum(u) / 2 .- U/2 # hide
FT_dat_pad = abs.(KomaMRI.fftc(dat_pad; dims=1)) # hide
scale_factor = maximum(abs.(M.xy)) / maximum(FT_dat_pad) # hide
s3 = scatter(x=u, y=FT_dat_pad*scale_factor, name="|FT(B₁(t))|", line=attr(dash="dash")) # hide
# Plot all results
pb = plot([s1,s2,s3], Layout(title="30 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]", xaxis_range=[-fmax,fmax])) # hide
if isFileMD savefig(pb, "../assets/4b-profile.html") end # hide
if isFileJL display(pb) end # hide
nothing # hide

#md # Which produces the follwoing slice profile:

#md # ```@raw html
#md # <object type="text/html" data="../../assets/4b-profile.html" style="width:100%; height:380px;"></object>
#md # ```

#md # For a flip angle of 30 deg the slice profile is very close to
#md # the small tip angle approximation, the Fourier transform of ``B_{1}(t)``.

#md # But what will happen if we use a flip angle of 120 deg instead?

α_desired = 120 + 0im               # The multiplication of a complex number scales the RF pulse of a Sequence
α = KomaMRI.get_flip_angles(seq)[1] #Previous FA approx 30 deg
seq = (α_desired / α) * seq         # Scaling the pulse to have a flip angle of 120
M = simulate_slice_profile(seq; z, simParams)

#Plot simulation results
s1 = scatter(x=f, y=abs.(M.xy), name="|Mxy|") # hide
# Small tip approximation
dat = seq.RF.A[1] # hide
N = length(dat) # hide
dat_pad = [zeros(floor(Int64,N)); dat; zeros(floor(Int64,N))] # hide
N_pad = length(dat_pad) # hide
U = 1 / (Trf) * N / N_pad #hide
u = range(0, (N_pad - 1) * U; step=U) # hide
u = u .- maximum(u) / 2 .- U/2 # hide
FT_dat_pad = abs.(KomaMRI.fftc(dat_pad; dims=1)) # hide
scale_factor = maximum(abs.(M.xy)) / maximum(FT_dat_pad) # hide
s2 = scatter(x=u, y=FT_dat_pad*scale_factor, name="|FT(B₁(t))|", line=attr(dash="dash")) # hide
# Plot all results
pa = plot([s1,s2], Layout(title="120 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]", xaxis_range=[-fmax,fmax])) # hide
if isFileMD savefig(pa, "../assets/4a-profile.html") end # hide
if isFileJL display(pa) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/4a-profile.html" style="width:100%; height:320px;"></object>
#md # ```

#md # For this case the small tip angle approximation breaks 😢 (thus, the reason for its name!), 
#md # and the slice profile is not similar to ``\mathcal{F}\left(B_{1}\left(t\right)\right)``. 

#md # This basic sinc pulse is not very ``B_{1}`` insensitive, for this, is better to use adiabatic RF pulses.
#md # Watch out for a future example 👀 with adiabatic RF pulses.
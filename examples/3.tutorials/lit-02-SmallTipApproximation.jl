# # Small Tip Angle Approximation

# > Based on the results in page 41 of the book "Handbook of MRI Pulse Sequences" by Bernstein et al.

using KomaMRI # hide
sys = Scanner() # hide
sys.limits.Smax = 50 # hide

# In this example, we will showcase a common approximation in MRI, the small tip angle approximation.
# For this, we will simulate a slice profile for spins with positions ``z\in[-2,\,2]\,\mathrm{cm}``
# and with a gradient ``G_{z}`` so their frequencies are mapped to ``f\in[-5,\,5]\,\mathrm{kHz}``.
# To start, we define an RF pulse with a flip angle of 30 deg and pulse duration of ``T_{\mathrm{rf}}=3.2\,\mathrm{ms}``.

B1 = 4.92e-6
Trf = 3.2e-3
zmax = 2e-2
fmax = 5e3
z = range(-zmax, zmax, 400)
Gz = fmax / (γ * zmax)
f = γ * Gz * z # hide

# The designed RF pulse is presented in the figure below,
# where the additional gradient refocuses the spins' phase after the excitation.

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
p2 = plot_seq(seq; max_rf_samples=Inf, slider=false)
#md savefig(p2, "../assets/42-seq.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/42-seq.html" style="width:100%; height:380px;"></object>
#md # ```

# Now we will perform the simulation using the function [`simulate_slice_profile`](@ref).
# Note that we modified `Δt_rf` in `sim_params` to match the resolution of the waveform.

sim_params = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, sim_params)

using PlotlyJS # hide
s1 = scatter(x=f, y=real.(M.xy), name="Mx") # hide
s2 = scatter(x=f, y=imag.(M.xy), name="My") # hide
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
pb = plot([s1,s2,s3], Layout(title="30 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]", xaxis_range=[-fmax,fmax])) # hide
#md savefig(pb, "../assets/4b-profile.html") # hide
#jl display(pb)

#md # This produces the following slice profile:

#md # ```@raw html
#md # <object type="text/html" data="../../assets/4b-profile.html" style="width:100%; height:380px;"></object>
#md # ```

# As you can see, for a flip angle of 30 deg, the slice profile is very close to
# the small tip angle approximation (the Fourier transform of ``B_{1}(t)``).

# But what will happen if we use a flip angle of 120 deg instead?

α_desired = 120 + 0im               # The multiplication of a complex number scales the RF pulse of a Sequence
α = get_flip_angles(seq)[1] # Previous FA approx 30 deg
seq = (α_desired / α) * seq         # Scaling the pulse to have a flip angle of 120
M = simulate_slice_profile(seq; z, sim_params)

s1 = scatter(x=f, y=abs.(M.xy), name="|Mxy|") # hide
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
pa = plot([s1,s2], Layout(title="120 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]", xaxis_range=[-fmax,fmax])) # hide
#md savefig(pa, "../assets/4a-profile.html") # hide
#jl display(pa)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/4a-profile.html" style="width:100%; height:320px;"></object>
#md # ```

# For this case, the small tip angle approximation breaks 😢, thus, the reason for its name!

# This basic sinc pulse is not designed to be ``B_{1}``-insensitive.  Some adiabatic RF pulses have been proposed to achieve this.
# Watch out for a future example showing these adiabatic RF pulses 👀.

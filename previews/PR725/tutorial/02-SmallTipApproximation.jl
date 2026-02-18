using KomaMRI, Suppressor #hide
sys = Scanner() #hide
sys.Smax = 50; #hide

B1 = 4.92e-6
Trf = 3.2e-3
zmax = 2e-2
fmax = 5e3
z = range(-zmax, zmax, 400)
Gz = fmax / (γ * zmax);
f = γ * Gz * z; #hide

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
p2 = plot_seq(seq; max_rf_samples=Inf, slider=false)
display(p2)

sim_params = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = @suppress simulate_slice_profile(seq; z, sim_params)

using PlotlyJS #hide
s1 = scatter(x=f, y=real.(M.xy), name="Mx") #hide
s2 = scatter(x=f, y=imag.(M.xy), name="My") #hide
dat = seq.RF.A[1] #hide
N = length(dat) #hide
dat_pad = [zeros(floor(Int64,N)); dat; zeros(floor(Int64,N))] #hide
N_pad = length(dat_pad) #hide
U = 1 / (Trf) * N / N_pad #hide
u = range(0, (N_pad - 1) * U; step=U) #hide
u = u .- maximum(u) / 2 .- U/2 #hide
FT_dat_pad = abs.(KomaMRI.fftc(dat_pad; dims=1)) #hide
scale_factor = maximum(abs.(M.xy)) / maximum(FT_dat_pad) #hide
s3 = scatter(x=u, y=FT_dat_pad*scale_factor, name="|FT(B₁(t))|", line=attr(dash="dash")) #hide
pb = plot([s1,s2,s3], Layout(title="30 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]", xaxis_range=[-fmax,fmax])) #hide
display(pb)

α_desired = 120 + 0im               # The multiplication of a complex number scales the RF pulse of a Sequence
α = get_flip_angles(seq)[1]         # Previous FA approx 30 deg
seq = (α_desired / α) * seq         # Scaling the pulse to have a flip angle of 120
M = @suppress simulate_slice_profile(seq; z, sim_params);

s1 = scatter(x=f, y=abs.(M.xy), name="|Mxy|") #hide
dat = seq.RF.A[1] #hide
N = length(dat) #hide
dat_pad = [zeros(floor(Int64,N)); dat; zeros(floor(Int64,N))] #hide
N_pad = length(dat_pad) #hide
U = 1 / (Trf) * N / N_pad #hide
u = range(0, (N_pad - 1) * U; step=U) #hide
u = u .- maximum(u) / 2 .- U/2 #hide
FT_dat_pad = abs.(KomaMRI.fftc(dat_pad; dims=1)) #hide
scale_factor = maximum(abs.(M.xy)) / maximum(FT_dat_pad) #hide
s2 = scatter(x=u, y=FT_dat_pad*scale_factor, name="|FT(B₁(t))|", line=attr(dash="dash")) #hide
pa = plot([s1,s2], Layout(title="120 deg SINC pulse (TBP=8, Hamming)", xaxis_title="Frequency [Hz]", xaxis_range=[-fmax,fmax])) #hide
display(pa)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

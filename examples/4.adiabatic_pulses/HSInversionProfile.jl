using KomaMRI, PlotlyJS
# RF Pulse Paramters, https://onlinelibrary.wiley.com/doi/epdf/10.1002/jmri.26021?saml_referrer
b1max = 13e-6      #Peak amplitude (uT)
Trf = 18.3e-3      #Pulse duration (ms)
β = 4e2            #frequency modulation param (rad/s)
μ = 6              #phase modulation parameter (dimensionless)
fmax = μ * β / (2π) # 2fmax = BW 
# Adiabatic condition b1max >> β*√μ/γ: 
b1max > β*sqrt(μ)/(2π*γ)
# Pulse Shape
t = range(-Trf/2, Trf/2, 201)
B1 =  b1max .* sech.(β.*t)
Δf = -fmax  .* tanh.(β.*t) 
# Sequence generation
fmax_sim = 2e3
Gz = fmax_sim / γ
seq = Sequence(
    [
        Grad(0.,0.);   #Gx
        Grad(0.,0.);   #Gy
        Grad(Gz,Trf,0) #Gz
    ;;], 
    [RF(B1,Trf,Δf,0);;]
    )
p1 = plot_seq(seq; max_rf_samples=Inf, slider=false)
KomaMRI.get_flip_angles(seq)[1]
# Simulation
simParams = Dict{String,Any}("Δt_rf"=>t[2]-t[1])
z = range(-1, 1, 400)
M = simulate_slice_profile(seq; simParams, z)
# Plot
f = γ*Gz*z
s1 = scatter(x=f,y=abs.(M.xy),name="|Mxy|")
s2 = scatter(x=f,y=M.z,name="Mz")
p2 = plot([s1,s2], Layout(title="Hyperbolic-Secant (HS) Adiabatic Inversion Pulse (μ=$μ, β=$β rad/s)", xaxis_title="Frequency [Hz]"))
[p1; p2]
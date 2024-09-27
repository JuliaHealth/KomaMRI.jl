# The purpose of this tutorial is to showcase the effect of the diffusion phenomenon over the acquisition
# of the MRI signal. In this script the simplest example is given, assuming isotropic diffusion, which
# can be characterized by only one coefficient. Theoretically, diffusion phenomenon measures the attenuation
# of the signal due to Brownian motion, a microscopical random movement.

# First necessary thing is a diffusion-sensitive sequence. In this case the selected sequence is one of the
# most typical, pulse gradient spin echo (PGSE). This sequence is characterized by the use of two diffusion
# gradients, placed right before and right after the inversion RF pulse. The duration of each gradient is
# defined by the δ parameter and the distance between the beginning of both gradients is described by the
# Δ parameter. In this tutorial a square shape will be assumed.

using KomaMRI
using LinearAlgebra
using Random

sys = Scanner()

durRF = π/2/(2π*γ*sys.B1)

rf = PulseDesigner.RF_hard(sys.B1,durRF,sys)
rf_inv = (0.0+2im)*rf

Δ = 43.1e-3
δ = 10.6e-3
τ = Δ - δ/3

G(b) = sqrt(b*1e6)/γ/2/π/δ/sqrt(τ)

# As can be seen, gradient amplitude is defined as a function of a 'b' parameter, which happens to be
# one of the most important parameters in order to describe the diffusion phenomenon, measuring the time
# per area (as opposed to diffusivity), given insight of the tissue at study. There is a 10^6 factor
# as this simulator use the International System of Units.

# Creation of the sequence

g_uni = Grad(1,δ)

seq = Sequence()
seq += rf
seq += Delay(2e-2)
seq += g_uni # Unitary gradient used in order to ease the sequence manipulation
seq += Delay((Δ-δ)/2)
seq += rf_inv
seq += Delay((Δ-δ)/2)
seq += g_uni # Unitary gradient used in order to ease the sequence manipulation
seq += Delay(2e-2)

adc = ADC(1,1e-6) # Analog/Digital converter (not apparent diffusion coefficient)

seq += adc

plot_seq(seq; show_adc = true) # Plotting the sequence 

# Now we have a "template" sequence that need to be adapted to meet the requirements of the desired acquisition
b0 = 0 # Defining a b-value for baseline acquisition
b = [300,600,900,1200,2000,3600] # Defining a b-value for diffusion acquisition
b_vector_x = [1,0,0] # Defining a b-vector to encode the direction of the gradient. (As is isotropic it makes no difference) must be ortonormal.
b_vector_y = [0,1,0] # Defining a b-vector to encode the direction of the gradient. (As is isotropic it makes no difference) must be ortonormal.
b_vector_z = [0,0,1] # Defining a b-vector to encode the direction of the gradient. (As is isotropic it makes no difference) must be ortonormal.
b_vector_eq = [sqrt(3)/3,sqrt(3)/3,sqrt(3)/3] # Defining a b-vector to encode the direction of the gradient. (As is isotropic it makes no difference) must be ortonormal.
R_x = hcat(b_vector_x,zeros(3,2)) # Vector used for rotating the gradient
R_y = hcat(b_vector_y,zeros(3,2)) # Vector used for rotating the gradient
R_z = hcat(b_vector_z,zeros(3,2)) # Vector used for rotating the gradient
R_eq = hcat(b_vector_eq,zeros(3,2)) # Vector used for rotating the gradient

seq_baseline = G(b0)*seq
seq_diffusion_x = []
seq_diffusion_y = []
seq_diffusion_z = []
seq_diffusion_eq = []
for i in b
    seq_x = R_x*G(i)*seq
    seq_y = R_y*G(i)*seq
    seq_z = R_z*G(i)*seq
    seq_eq = R_eq*G(i)*seq
    push!(seq_diffusion_x,seq_x)
    push!(seq_diffusion_y,seq_y)
    push!(seq_diffusion_z,seq_z)
    push!(seq_diffusion_eq,seq_eq)
end

plot_seq(seq_baseline)
plot_seq(seq_diffusion_x[1])
plot_seq(seq_diffusion_y[1])
plot_seq(seq_diffusion_z[1])
plot_seq(seq_diffusion_eq[1])

# With the sequence set, we have to define a phantom_iso. For this example, a collection of 100 spins grouped 
# in the same point will suffice the requirements. 

phantom_iso = Phantom(x=collect(range(0,0,1000)))
phantom_x= Phantom(x=collect(range(0,0,1000)))

# Typical values
phantom_iso.ρ .= 1
phantom_iso.T1 .= 1000e-3
phantom_iso.T2 .= 100e-3

adc_iso = sqrt(2e-9)*Diagonal([1.,1.,1.]) # Diffusivity coefficient in the main diagonal of a matrix
adc_x = sqrt(2e-9)*Diagonal([1.,0.,0.]) # Diffusivity in the x-component of the matrix

Nspins = length(phantom_iso)
Nt = 100
Δt = dur(seq)/Nt

# Initializing displacement directions 
dx = zeros(Nspins,Nt)
dy = zeros(Nspins,Nt)
dz = zeros(Nspins,Nt)
rng = MersenneTwister(1234) #Random seed
for i = 1:Nt
    tmp = sqrt(2*Δt)*adc_iso*randn(rng,3,Nspins) # Isotropic diffusivity
    dx[:,i] = tmp[1,:]
    dy[:,i] = tmp[2,:]
    dz[:,i] = tmp[3,:]
end
# Diffusion motion is the cumulative sum of displacements
dx = cumsum(dx;dims=2)
dy = cumsum(dy;dims=2)
dz = cumsum(dz;dims=2)
phantom_iso.motion = MotionList(Path(dx,dy,dz,TimeRange(0.,dur(seq)))) # Adding the Brownian motion to the spins
plot_phantom_map(phantom_iso,:T1;intermediate_time_samples=100) # If used in Visual Studio, should be plotted outside of it

# Initializing displacement directions 
dx = zeros(Nspins,Nt)
dy = zeros(Nspins,Nt)
dz = zeros(Nspins,Nt)
rng = MersenneTwister(1234) #Random seed
for i = 1:Nt
    tmp = sqrt(2*Δt)*adc_x*randn(rng,3,Nspins) # Diffusion only on the x-axis
    dx[:,i] = tmp[1,:]
    dy[:,i] = tmp[2,:]
    dz[:,i] = tmp[3,:]
end
# Diffusion motion is the cumulative sum of displacements
dx = cumsum(dx;dims=2)
dy = cumsum(dy;dims=2)
dz = cumsum(dz;dims=2)
phantom_x.motion = MotionList(Path(dx,dy,dz,TimeRange(0.,dur(seq)))) # Adding the Brownian motion to the spins
plot_phantom_map(phantom_x,:T1;intermediate_time_samples=100) # If used in Visual Studio, should be plotted outside of it
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"

# Simulating diffusion sequences
S0_iso = simulate(phantom_iso, seq_baseline, sys; sim_params) # Simulation of the baseline acquisition
S1_iso_x = ComplexF32[]
S1_iso_y = ComplexF32[]
S1_iso_z = ComplexF32[]
S1_iso_eq = ComplexF32[]
for i = 1:length(b)
    S1_x = simulate(phantom_iso, seq_diffusion_x[i], sys; sim_params) # Simulation of the diffusion acquisition
    S1_y = simulate(phantom_iso, seq_diffusion_y[i], sys; sim_params) # Simulation of the diffusion acquisition
    S1_z = simulate(phantom_iso, seq_diffusion_z[i], sys; sim_params) # Simulation of the diffusion acquisition
    S1_eq = simulate(phantom_iso, seq_diffusion_eq[i], sys; sim_params) # Simulation of the diffusion acquisition
    push!(S1_iso_x,S1_x[1,1,1])
    push!(S1_iso_y,S1_y[1,1,1])
    push!(S1_iso_z,S1_z[1,1,1])
    push!(S1_iso_eq,S1_eq[1,1,1])
end
E_iso_x = abs.(S1_iso_x)./abs.(S0_iso[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline
E_iso_y = abs.(S1_iso_y)./abs.(S0_iso[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline
E_iso_z = abs.(S1_iso_z)./abs.(S0_iso[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline
E_iso_eq = abs.(S1_iso_eq)./abs.(S0_iso[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline

E_iso_theoretical = exp.(-b.*1e6.*2e-9)

S0_x = simulate(phantom_x, seq_baseline, sys; sim_params) # Simulation of the baseline acquisition
S1_x_x = ComplexF32[]
S1_x_y = ComplexF32[]
S1_x_z = ComplexF32[]
S1_x_eq = ComplexF32[]
for i = 1:length(b)
    S1_x = simulate(phantom_x, seq_diffusion_x[i], sys; sim_params) # Simulation of the diffusion acquisition
    S1_y = simulate(phantom_x, seq_diffusion_y[i], sys; sim_params) # Simulation of the diffusion acquisition
    S1_z = simulate(phantom_x, seq_diffusion_z[i], sys; sim_params) # Simulation of the diffusion acquisition
    S1_eq = simulate(phantom_x, seq_diffusion_eq[i], sys; sim_params) # Simulation of the diffusion acquisition
    push!(S1_x_x,S1_x[1,1,1])
    push!(S1_x_y,S1_y[1,1,1])
    push!(S1_x_z,S1_z[1,1,1])
    push!(S1_x_eq,S1_eq[1,1,1])
end
E_x_x = abs.(S1_x_x)./abs.(S0_x[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline
E_x_y = abs.(S1_x_y)./abs.(S0_x[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline
E_x_z = abs.(S1_x_z)./abs.(S0_x[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline
E_x_eq = abs.(S1_x_eq)./abs.(S0_x[1,1,1]) # Computing attenuation as the ratio between diffusion and baseline

## Plotting isotropic attenuation

using Plots

plot([300,600,900,1200,2000,3600],E_iso_x,label="Attenuation with [1,0,0] b-vector",title="Isotropic diffusivity")
plot!([300,600,900,1200,2000,3600],E_iso_y,label="Attenuation with [0,1,0] b-vector")
plot!([300,600,900,1200,2000,3600],E_iso_z,label="Attenuation with [0,0,1] b-vector")
plot!([300,600,900,1200,2000,3600],E_iso_eq,label="Attenuation with [√3/3,√3/3,√3/3] b-vector")
plot!([300,600,900,1200,2000,3600],E_iso_theoretical,label="Theoretical attenuation")

# In this figure all sequences get similar results since the placement of the gradients does not affect
# the received signal. The random nature of Brownian motion makes the differences in signal loss in the different
# axis.

## Plotting x-axis attenuation

plot([300,600,900,1200,2000,3600],E_x_x,label="Attenuation with [1,0,0] b-vector",title="X-axis diffusivity",xlabel = "b-value",ylabel = "attenuation")
plot!([300,600,900,1200,2000,3600],E_x_y,label="Attenuation with [0,1,0] b-vector")
plot!([300,600,900,1200,2000,3600],E_x_z,label="Attenuation with [0,0,1] b-vector")
plot!([300,600,900,1200,2000,3600],E_x_eq,label="Attenuation with [√3/3,√3/3,√3/3] b-vector")
plot!([300,600,900,1200,2000,3600],E_iso_theoretical,label="Theoretical attenuation")

# In this figure, y and z b-vectors do not get diffusion signal at all, as they remain in 1. 
# The gradient placed in the x axis receive an amount of signal similar to the first scenario, relatable to the theoretical value.
# The gradient placed in the equidistant point suffers much less attenuation.

# On both figures the effect of the b-value is patent, exponentially reducing the signal as it gets bigger.
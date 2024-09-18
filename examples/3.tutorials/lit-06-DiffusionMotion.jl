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

Δ = .0735
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
b = 300 # Defining a b-value for diffusion acquisition
b_vector = [1,0,0] # Defining a b-vector to encode the direction of the gradient. (As is isotropic it makes no difference) must be ortonormal.
R = hcat(b_vector,zeros(3,2)) # Vector used for rotating the gradient

seq_baseline = G(b0)*seq
seq_diffusion = R*G(b)*seq

plot_seq(seq_baseline)
plot_seq(seq_diffusion)

# With the sequence set, we have to define a phantom. For this example, a collection of 100 spins grouped 
# in the same point will suffice the requirements. 

phantom = Phantom(x=collect(range(0,0,1000)))

# Typical values
phantom.ρ .= 1
phantom.T1 .= 1000e-3
phantom.T2 .= 100e-3

adc_iso = sqrt(2e-9)*Diagonal([1.,1.,1.]) # Diffusivity coefficient in the main diagonal of a matrix

Nspins = length(phantom)
Nt = 100
Δt = dur(seq)/Nt

# Initializing displacement directions 
dx = zeros(Nspins,Nt)
dy = zeros(Nspins,Nt)
dz = zeros(Nspins,Nt)
rng = MersenneTwister(1234) #Random seed
for i = 1:Nt
    tmp = sqrt(2*Δt)*adc_iso*randn(rng,3,Nspins) # Assuming isotropic gaussian motion
    dx[:,i] = tmp[1,:]
    dy[:,i] = tmp[2,:]
    dz[:,i] = tmp[3,:]
end
# Diffusion motion is the cumulative sum of displacements
dx = cumsum(dx;dims=2)
dy = cumsum(dy;dims=2)
dz = cumsum(dz;dims=2)
phantom.motion = ArbitraryMotion(0.,dur(seq), dx, dy, dz) # Adding the Brownian motion to the spins
plot_phantom_map(phantom,:T1) # If used in Visual Studio, should be plotted outside of it

sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
S0 = simulate(phantom, seq_baseline, sys; sim_params) # Simulation of the baseline acquisition
S1 = simulate(phantom, seq_diffusion, sys; sim_params) # Simulation of the diffusion acquisition
E = abs.(S1)./abs.(S0) # Computing attenuation as the ratio between diffusion and baseline

E_theoric = exp(-b*1e6*2e-9) # Checking the result

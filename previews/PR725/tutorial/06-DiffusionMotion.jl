using KomaMRI #hide
using PlotlyJS #hide
using Random, Suppressor #hide

Nspins = 10_000
obj = Phantom(;
    x  = zeros(Nspins),
    T1 = ones(Nspins) * 1000e-3,
    T2 = ones(Nspins) * 100e-3,
);

D = 2e-9                # Diffusion Coefficient of water in m^2/s
T = 100e-3              # Duration of the motion
Nt = 100                # Number of time steps
Δt = T / (Nt - 1)       # Time sep
Δr = sqrt(2 * D * Δt);  # √ Mean square displacement

rng = MersenneTwister(1234) # Setting up the random seed
dx = cumsum([zeros(Nspins) Δr .* randn(rng, Nspins, Nt - 1)]; dims=2)
dy = cumsum([zeros(Nspins) Δr .* randn(rng, Nspins, Nt - 1)]; dims=2)
dz = cumsum([zeros(Nspins) Δr .* randn(rng, Nspins, Nt - 1)]; dims=2);

random_walk = KomaMRI.path(dx, dy, dz, TimeRange(0.0, T))
obj.motion = random_walk
p1 = plot_phantom_map(obj, :T1; time_samples=Nt÷4, height=450)
display(p1)

sys   = Scanner()
durRF = 1e-3
B1    = (π / 2) / (2π * γ * durRF)
rf90  = PulseDesigner.RF_hard(B1, durRF, sys)
rf180 = (0.0 + 2im) * rf90;

G = 30e-3            # Gradient amplitude
δ = 30e-3            # Duration of the gradient
Δ = durRF + δ        # Time between the two gradients
gx_diff = Grad(G, δ);

adc_dwell_time = 1e-6
adc = ADC(1, adc_dwell_time, durRF/2 - adc_dwell_time/2); # ADCs with N=1 are positioned at the center

seq = Sequence()
seq += rf90
seq += gx_diff
seq += rf180
seq += gx_diff
seq += adc
p2 = plot_seq(seq; show_adc=true) # Plotting the sequence
display(p2);

function bvalue(seq)
    block, axis = 2, 1 # Gx from second block
    G = seq.GR[axis, block].A
    δ = seq.GR[axis, block].T
    Δ = dur(seq[2:3]) # Because there are no gaps
    b = (2π * γ * G * δ)^2 * (Δ - δ/3)
    return b * 1e-6
end;

seqs = Sequence[] # Vector of sequences
bvals = [0, 250, 500, 1000, 1500, 2000] # b-values in s/mm^2
for bval_target in bvals
    gradient_scaling = sqrt(bval_target / bvalue(seq))
    seq_b = gradient_scaling * seq
    push!(seqs, seq_b)
end

sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["Δt"] = Δt; # Set max. grad. time step to fit diffusion time step

signals = @suppress simulate.(Ref(obj), seqs, Ref(sys); sim_params); # simulate broadcasted over seqs

Sb = [sb[1] for sb in signals] # Reshaping the simulated signals
bvals_si = bvals .* 1e6; # Convert b-values from s/mm^2 to s/m^2

E_simulated   = abs.(Sb) ./ abs.(Sb[1])
E_theoretical = exp.(-bvals_si .* D);
s_sim  = scatter(x=bvals, y=E_simulated,   name="Simulated") #hide
s_theo = scatter(x=bvals, y=E_theoretical, name="exp(-b D)", line=attr(dash="dash")) #hide
layout = Layout(title="Diffusion-induced signal attenuation E(b)", xaxis=attr(title="b-value [s/mm^2]")) #hide
p3 = plot([s_sim, s_theo], layout); #hide
display(p3);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

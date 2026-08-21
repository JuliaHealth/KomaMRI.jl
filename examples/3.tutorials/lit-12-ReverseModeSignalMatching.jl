## Estimating Spin Density from an FID

using Enzyme
using KomaMRI
using Reactant

Reactant.set_default_backend("cpu")

# After a hard ``90^\circ`` pulse, spins with different off-resonance frequencies contribute distinct phase evolutions to the measured free induction decay (FID). As a result, we can identify their spin densities from the complex signal.
#
# Let ``A(\rho)`` be the simulated FID for spin densities ``\rho`` and let
# ``b`` be the measured signal. We recover the densities by minimizing
#
# ```math
# L(\rho) = \frac{1}{N}\left\|A(\rho)-b\right\|_2^2.
# ```

## Defining the FID experiment
#
# A hard ``90^\circ`` pulse creates transverse magnetization, and we will samples its FID at nine points.

Trf = 0.4e-3
B1 = (π / 2) / (2π * γ * Trf)
seq = Sequence()
@addblock seq += RF(B1, Trf)
@addblock seq += ADC(9, 4e-3)

# We'll use three spins with distinct off-resonance frequencies, so each spin will contribute a different complex oscillation to our measured signal.

obj = Phantom(;
    x=zeros(3),
    ρ=ones(3),
    T1=ones(3),
    T2=fill(80e-3, 3),
    Δw=2π .* [-120.0, 0.0, 120.0],
)
sim_params = Dict{String,Any}(
    "gpu" => false,
    "Nthreads" => 1,
    "return_type" => "mat",
    "Δt" => 50e-6,
    "Δt_rf" => 25e-6,
)
params = (; seq, obj, sys=Scanner(), sim_params)

## Defining the forward model and synthetic measurement
#
# `forward` replaces the three spin densities and simulates the FID.

function forward(ρ, p)
    obj = copy(p.obj)
    obj.ρ .= ρ
    return simulate(obj, p.seq, p.sys; sim_params=p.sim_params, verbose=false)
end

# In this case we're generating a target ('b') based on spin density values we want the optimization to recover

b = forward([0.3, 0.8, 0.5], params)
params = merge(params, (; b))

loss(ρ, p) = sum(abs2, forward(ρ, p) .- p.b) / length(p.b)

## Differentiating the signal model
#
# Enzyme differentiates through `simulate` and returns the sensitivity of the signal mismatch to each spin density.

function loss_and_gradient(ρ, p)
    result = Enzyme.gradient(
        Enzyme.ReverseWithPrimal,
        loss,
        ρ,
        Enzyme.Const(p),
    )
    return result.val, result.derivs[1]
end

## Recovering the spin densities
#
# Reactant compiles the forward simulation and Enzyme reverse pass as a single reusable function. In this case, only the phantom and target signal need to be traced.

params_ra = merge(params, (;
    obj=Reactant.to_rarray(params.obj),
    b=Reactant.to_rarray(params.b),
))
ρ = Reactant.to_rarray(fill(0.7, 3))
compiled = Reactant.@allowscalar Reactant.compile(
    loss_and_gradient,
    (ρ, params_ra);
    sync=true,
)

# We'll use 20 gradient-descent steps to move from the uniform initial guess toward the proton densities used to generate `b`.

initial_loss = Reactant.to_number(first(compiled(ρ, params_ra)))
for _ in 1:20
    _, ∇loss = compiled(ρ, params_ra)
    global ρ = ρ .- 0.15 .* ∇loss
end
final_loss = Reactant.to_number(first(compiled(ρ, params_ra)))

(;
    initial_loss=round(initial_loss; sigdigits=4),
    final_loss=round(final_loss; sigdigits=4),
    recovered_density=round.(Array(ρ); digits=3),
)

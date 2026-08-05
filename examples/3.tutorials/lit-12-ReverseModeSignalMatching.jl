# # Reverse-mode signal matching
#
# The same Reactant and Enzyme path can optimize parameters that affect an
# acquired ADC signal. In this example, reverse-mode AD recovers three spin
# densities from a short complex FID.

using Enzyme
using KomaMRI
using Reactant

Reactant.set_default_backend("cpu")
Reactant.allowscalar(false)

# A hard 90-degree pulse is followed by nine ADC samples. Distinct
# off-resonance frequencies make each spin's contribution identifiable.
rf_duration = 0.4e-3
b1 = (π / 2) / (2π * γ * rf_duration)
seq = Sequence()
@addblock seq += RF(b1, rf_duration)
@addblock seq += ADC(9, 4e-3)

obj = Phantom(;
    x=zeros(3),
    ρ=ones(3),
    T1=ones(3),
    T2=fill(80e-3, 3),
    Δw=2π .* [-120.0, 0.0, 120.0],
)
sim_params = Dict{String,Any}(
    "sim_method" => Bloch(),
    "gpu" => false,
    "Nthreads" => 1,
    "return_type" => "mat",
    "precision" => "f64",
    "sampling_rule" => MaxStepSizeRule(50e-6, 25e-6),
)
params = (; seq, obj, sys=Scanner(), sim_params)

function signal_forward(density, params)
    obj_aux = copy(params.obj)
    obj_aux.ρ .= density
    return simulate(
        obj_aux,
        params.seq,
        params.sys;
        sim_params=params.sim_params,
        verbose=false,
    )
end

target_density = [0.3, 0.8, 0.5]
target_signal = signal_forward(target_density, params)
params = merge(params, (; target_signal))

function signal_loss(density, params)
    signal = signal_forward(density, params)
    return sum(abs2, signal .- params.target_signal) / length(signal)
end

function signal_loss_and_gradient(density, params)
    result = Enzyme.gradient(
        Enzyme.ReverseWithPrimal,
        signal_loss,
        density,
        Enzyme.Const(params),
    )
    return result.val, result.derivs[1]
end

params_device = merge(params, (;
    obj=Reactant.to_rarray(params.obj),
    target_signal=Reactant.to_rarray(params.target_signal),
))
density = Reactant.to_rarray(fill(0.7, 3))
compiled_signal_gradient = Reactant.@compile sync=true signal_loss_and_gradient(
    density,
    params_device,
)

function optimize_signal(density, params, gradient_function; iterations=12, step_size=0.15)
    initial_loss = Reactant.to_number(first(gradient_function(density, params)))
    for _ in 1:iterations
        _, gradient = gradient_function(density, params)
        density = density .- step_size .* gradient
    end
    final_loss = Reactant.to_number(first(gradient_function(density, params)))
    return density, initial_loss, final_loss
end

density, initial_loss, final_loss = optimize_signal(
    density,
    params_device,
    compiled_signal_gradient,
)

(;
    initial_loss=round(initial_loss; sigdigits=4),
    final_loss=round(final_loss; sigdigits=4),
    recovered_density=round.(Array(density); digits=3),
)

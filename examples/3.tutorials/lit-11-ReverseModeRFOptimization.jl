# # Reverse-mode RF pulse optimization
#
# Optimize three RF control points through interpolation and a Bloch simulation
# using Enzyme reverse mode compiled by Reactant.

using Enzyme
using KomaMRI
using Reactant

Reactant.set_default_backend("cpu")
Reactant.allowscalar(false)

Trf = 0.6e-3
seq = Sequence()
@addblock seq += (RF(zeros(ComplexF64, 7), Trf), z=Grad(8e-3, Trf))

params = (;
    seq,
    obj=Phantom(; x=zeros(9), z=collect(range(-4e-3, 4e-3; length=9))),
    sys=Scanner(),
    sim_params=Dict{String,Any}(
        "sim_method" => Bloch(),
        "gpu" => false,
        "Nthreads" => 1,
        "return_type" => "state",
        "precision" => "f64",
        "sampling_rule" => MaxStepSizeRule(50e-6, 25e-6),
    ),
    node_times=range(0.0, Trf; length=3),
    sample_times=range(0.0, Trf; length=7),
)

# The optimization variables are RF nodes in microtesla. They are interpolated
# to the samples used by `simulate` inside the differentiated forward model.
function rf_forward(nodes, p)
    samples = KomaMRIBase.linear_interpolate_samples(
        (t=p.node_times, A=nodes),
        p.sample_times,
    )
    seq = KomaMRIBase.set_rf_amplitude(p.seq, 1e-6 .* complex.(samples))
    return simulate(p.obj, seq, p.sys; sim_params=p.sim_params, verbose=false).xy
end

target = rf_forward([0.0, 5.0, 0.0], params)
params = merge(params, (; target))

rf_loss(nodes, p) = sum(abs2, rf_forward(nodes, p) .- p.target) / length(p.target)

function rf_loss_and_gradient(nodes, p)
    result = Enzyme.gradient(
        Enzyme.ReverseWithPrimal,
        rf_loss,
        nodes,
        Enzyme.Const(p),
    )
    return result.val, result.derivs[1]
end

# Transfer only the arrays used by the compiled simulation.
params_device = merge(params, (;
    obj=Reactant.to_rarray(params.obj),
    target=Reactant.to_rarray(params.target),
))
nodes = Reactant.to_rarray([0.0, 3.5, 0.0])
compiled_gradient = Reactant.@compile sync=true rf_loss_and_gradient(nodes, params_device)

initial_loss = Reactant.to_number(first(compiled_gradient(nodes, params_device)))
for _ in 1:12
    _, ∇loss = compiled_gradient(nodes, params_device)
    global nodes = nodes .- 8.0 .* ∇loss
end
final_loss = Reactant.to_number(first(compiled_gradient(nodes, params_device)))

(;
    initial_loss=round(initial_loss; sigdigits=4),
    final_loss=round(final_loss; sigdigits=4),
    optimized_nodes=round.(Array(nodes); digits=3),
)

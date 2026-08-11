# # Designing a Slice-selective Excitation Profile

using Enzyme
using KomaMRI
using Reactant

Reactant.set_default_backend("cpu")

# A slice-selective gradient maps position to resonance frequency, so the RF waveform determines the transverse magnetization profile after excitation.
# Here, ``x`` contains the value of three RF control points in microtesla.
# These control points are interpolated to a finer RF timeline before simulation.
#
# If ``A(x)`` denotes the Bloch simulation and ``b`` is the desired profile, we
# design the optimal pulse by minimizing the objective function:
#
# ```math
# L(x) = \frac{1}{N}\left\|A(x)-b\right\|_2^2.
# ```

# ## Defining the slice-selection experiment
#
# The sequence contains a short RF pulse played with a constant slice-selection gradient.

Trf = 0.6e-3
seq = Sequence()
@addblock seq += (RF(zeros(ComplexF64, 7), Trf), z=Grad(8e-3, Trf))

# We place 9 spins along ``z`` to let us evaluate the excited slice profile.

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
        "Δt" => 50e-6,
        "Δt_rf" => 25e-6,
    ),
    node_times=range(0.0, Trf; length=3),
    sample_times=range(0.0, Trf; length=7),
)

# ## Defining the forward model and target
#
# `forward` interpolates the three controls onto the seven RF samples and then simulates the resulting pulse. It returns the final transverse magnetization at each location.

function forward(x, p)
    samples = KomaMRIBase.linear_interpolate_samples(
        (t=p.node_times, A=x),
        p.sample_times,
    )
    seq = KomaMRIBase.set_rf_amplitude(p.seq, 1e-6 .* complex.(samples))
    return simulate(p.obj, seq, p.sys; sim_params=p.sim_params, verbose=false).xy
end

# In this example we're generating a target profile from a known three-node pulse. However, `b` could be any kind of desired or measured excitation profile.

b = forward([0.0, 5.0, 0.0], params)
params = merge(params, (; b))

loss(x, p) = sum(abs2, forward(x, p) .- p.b) / length(p.b)

# ## Differentiating the simulation
#
# Enzyme differentiates the loss function and returns the loss value and its gradient w.r.t the RF control points

function loss_and_gradient(x, p)
    result = Enzyme.gradient(
        Enzyme.ReverseWithPrimal,
        loss,
        x,
        Enzyme.Const(p),
    )
    return result.val, result.derivs[1]
end

# ## Optimizing the RF
#
# Reactant compiles the loss and its Enzyme reverse pass together so the same function can be reused at every optimization step.
params_ra = merge(params, (;
    obj=Reactant.to_rarray(params.obj),
    b=Reactant.to_rarray(params.b),
))
x = Reactant.to_rarray([0.0, 3.5, 0.0])
compiled = Reactant.@allowscalar Reactant.compile(
    loss_and_gradient,
    (x, params_ra);
    sync=true,
)

# In this example we take 12 gradient descent steps toward a pulse that achieves the target

initial_loss = Reactant.to_number(first(compiled(x, params_ra)))
for _ in 1:12
    _, ∇loss = compiled(x, params_ra)
    global x = x .- 8.0 .* ∇loss
end
final_loss = Reactant.to_number(first(compiled(x, params_ra)))

(;
    initial_loss=round(initial_loss; sigdigits=4),
    final_loss=round(final_loss; sigdigits=4),
    optimized_nodes=round.(Array(x); digits=3),
)

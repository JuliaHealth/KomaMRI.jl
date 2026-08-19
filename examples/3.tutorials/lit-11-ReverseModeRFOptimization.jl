# # Designing a Slice-selective Excitation Profile

using Enzyme
using KomaMRI
using Reactant

Reactant.set_default_backend("cpu")

# A slice-selective gradient maps position to resonance frequency, so the RF waveform determines the transverse magnetization profile after excitation.
# Here, ``x`` contains the value of seven RF control points in microtesla.
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
# The sequence contains a short RF pulse played with a constant slice-selection gradient, followed by a gradient lobe that refocuses the transverse phase.

Trf = 2.4e-3
Gz = 8e-3
seq = Sequence()
@addblock seq += (RF(zeros(ComplexF64, 25), Trf), z=Grad(Gz, Trf))
@addblock seq += (z=Grad(-Gz, Trf / 2),)

# We place spins along ``z`` and prescribe a one-dimensional top-hat target with a 2 mm slice thickness.

z = collect(range(-4e-3, 4e-3; length=41))
b = @. im * Float64(abs(z) <= 1e-3)
params = (;
    seq,
    obj=Phantom(; x=zeros(length(z)), z),
    sys=Scanner(),
    sim_params=Dict{String,Any}(
        "gpu" => false,
        "Nthreads" => 1,
        "return_type" => "state",
        "Δt" => 100e-6,
        "Δt_rf" => 100e-6,
    ),
    node_times=range(0.0, Trf; length=7),
    sample_times=range(0.0, Trf; length=25),
    b,
)

# ## Defining the forward model and target
#
# `forward` interpolates the seven controls onto the RF waveform and then simulates the resulting pulse. It returns the final transverse magnetization at each location.

function forward(x, p)
    samples = KomaMRIBase.linear_interpolate_samples(
        (t=p.node_times, A=x),
        p.sample_times,
    )
    seq = copy(p.seq)
    seq.RF[1].A = 1e-6 .* complex.(samples)
    return simulate(p.obj, seq, p.sys; sim_params=p.sim_params, verbose=false).xy
end

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
seq_ra = Sequence()
@addblock seq_ra += (Reactant.to_rarray(params.seq.RF[1]), z=params.seq.GR[3, 1])
@addblock seq_ra += (z=params.seq.GR[3, 2],)
params_ra = merge(params, (;
    seq=seq_ra,
    obj=Reactant.to_rarray(params.obj),
    b=Reactant.to_rarray(params.b),
))
x = Reactant.to_rarray(zeros(7))
compiled = Reactant.@allowscalar Reactant.compile(
    loss_and_gradient,
    (x, params_ra);
    sync=true,
)

# In this example we take 20 gradient descent steps toward a pulse that achieves the target.

initial_loss = Reactant.to_number(first(compiled(x, params_ra)))
for _ in 1:20
    _, ∇loss = compiled(x, params_ra)
    global x = x .- 8.0 .* ∇loss
end
final_loss = Reactant.to_number(first(compiled(x, params_ra)))

# The optimized transverse magnetization can now be compared directly with the desired slice profile.

optimized_profile = forward(Array(x), params)
using PlotlyBase #hide
target_trace = scatter(x=z .* 1e3, y=abs.(b), name="Target", line=attr(dash="dash")) #hide
optimized_trace = scatter(x=z .* 1e3, y=abs.(optimized_profile), name="Optimized") #hide
profile_plot = Plot([target_trace, optimized_trace], Layout(xaxis_title="z [mm]", yaxis_title="|Mxy|", title="Designed slice profile")) #hide
#jl display(profile_plot)
#md profile_plot #hide

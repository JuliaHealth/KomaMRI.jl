# # Designing an Optimal 90° Excitation Pulse

# Designing a slice-selective RF pulse becomes difficult when include effects in our model that are typically assumed away (nonlinear Bloch dynamics, off-resonance, relaxation).
#
# In this tutorial, KomaMRI will be used for RF pulse design, showcasing its
# differentiable capabilities. We use Enzyme to differentiate a complete
# KomaMRI simulation and Reactant to compile the resulting gradient calculation.
# The objective is an optimized 90° excitation with a Butterworth slice profile.

using Enzyme, Reactant
using KomaMRIBase: cubic_interpolate_samples
Reactant.set_default_backend("cpu") # Change to "gpu" to compile for a GPU.
using KomaMRI #hide
using KomaMRI.PulseDesigner #hide
using PlotlyBase #hide
using Printf #hide

# ## Setting up the sequence
#
# We begin with a common 90° sinc pulse and its slice-selection and rephasing
# gradients. The RF lasts 3 ms, excites a 2 mm slice, and uses a 30 mT/m
# slice-selection gradient. Optimization will replace its waveform without
# rebuilding the sequence.

Trf = 3e-3
Gz = 30e-3
slice_thickness = 2e-3
sys = Scanner()
seq = build_sinc_pulse(
    π / 2;
    duration=Trf,
    slice_thickness,
    time_bw_product=γ * Gz * slice_thickness * Trf,
    apodization=0.46,
    sys,
    use=Excitation(),
)
## Reuse the sinc RF sampling grid for the optimized waveform.
rf = seq.RF[1];

# Optimizing every RF sample independently would be unnecessary. Instead, the
# design vector ``\boldsymbol{x}=(x_1,\ldots,x_P)^\mathsf{T}`` contains ``P=30``
# complex RF amplitudes in microtesla. The complex phase is part of
# ``\boldsymbol{x}`` itself. Every node is initialized to complex zero.

num_nodes = 30
sample_times = range(rf.delay, dur(rf); length=length(rf.A))
node_times = range(rf.delay, dur(rf); length=num_nodes);

# ## [Target magnetization profile ``\boldsymbol{b}``](@id target-magnetization-profile)
#
# The optimization grid contains ``N_{\mathrm{spins}}=41`` spins at locations
# ``z_j``. We define the desired complex transverse magnetization directly as
# the discrete vector
#
# ```math
# \boldsymbol{b}
# = \mathrm{i}\left[
#     \frac{1}{\sqrt{1 + (|z_j|/z_c)^{2n}}}
#   \right]_{j=1}^{N_{\mathrm{spins}}},
# ```
#
# with Butterworth degree ``n=8`` and cutoff ``z_c=1\,\mathrm{mm}``. Thus the
# target has unit ``M_y`` at the center and reaches half power at ``|z|=z_c``.

z_extent = 4e-3
Nspins = 41
butterworth_degree = 8
z_c = 1e-3 # Spatial cutoff at ±1 mm.
z = collect(range(-z_extent, z_extent; length=Nspins));
# Purely imaginary so the target lies entirely along My.
b = @. im / sqrt(1 + (abs(z) / z_c)^(2 * butterworth_degree));

# We use a deliberately coarse ``\Delta t_{\mathrm{rf}}=100\,\mu\mathrm{s}``
# during optimization so that differences between the simulation methods remain
# visible.

Δt_rf = 100e-6
p = (; # These parameters will be constant for the optimization.
    seq_template=seq,
    obj=Phantom(; x=zeros(Nspins), z),
    sys,
    sim_params=Dict{String,Any}(
        "Nthreads" => 1,
        "return_type" => "state",
        "Δt_rf" => Δt_rf,
    ),
    node_times,
    sample_times,
    b,
);

# ## [Forward operator ``\mathcal{A}_p(\boldsymbol{x})``](@id forward-operator)
#
# Let ``p`` collect everything held fixed during one optimization (the baseline
# sequence, spins, and simulation parameters). For a fixed
# ``p``, the forward calculation follows one path from the design variables to
# the profile that will be compared with ``\boldsymbol{b}``:
#
# ```math
# \boldsymbol{x}
# \xrightarrow{\text{interp. at }t_k}
# \left[B_1(t_k;\boldsymbol{x})\right]_{k=1}^{N_{\mathrm{rf}}}
# \xrightarrow[\mathcal{A}_p(\boldsymbol{x})]{\text{Simulation}}
# \boldsymbol{M}_{xy}.
# ```
#
# In code, `B1` contains the interpolated samples converted from microtesla to
# tesla and replaces only the amplitude vector of the first RF event. KomaMRI
# then returns one complex transverse-magnetization value for every spin on the
# design grid.

function 𝓐p(x, p)
    B1 = cubic_interpolate_samples((t=p.node_times, A=1e-6 .* x), p.sample_times)
    seq_x = copy(p.seq_template)
    seq_x.RF[1].A = B1
    Mxy = simulate(p.obj, seq_x, p.sys; sim_params=p.sim_params, verbose=false).xy
    return Mxy
end;

# ## [Loss function ``\mathcal{L}_p(\boldsymbol{x})``](@id loss-function)
#
# We minimize the mean squared error of the full complex profile:
#
# ```math
# \mathcal{L}_p(\boldsymbol{x})
# = \left\|\mathcal{A}_p(\boldsymbol{x})-\boldsymbol{b}\right\|_2^2 / N_{\mathrm{spins}}.
# ```
#
# Its Julia definition mirrors the equation:

𝓛p(x, p) = sum(abs2, 𝓐p(x, p) .- p.b) / length(p.b);

# ## [Differentiating ``\mathcal{L}_p(\boldsymbol{x})``](@id differentiating-loss)
#
# [Enzyme](https://enzyme.mit.edu/julia/stable/) applies reverse-mode automatic
# differentiation to the scalar loss and obtains all 30
# derivatives in one reverse pass. [Reactant](https://enzymead.github.io/Reactant.jl/stable/introduction/)
# traces this Julia calculation, lowers it through MLIR/XLA, performs targeted optimizations, and produces a
# compiled function that can be reused by every optimization step.

∇𝓛p(x, p) = Enzyme.gradient(Enzyme.Reverse, 𝓛p, x, Enzyme.Const(p))[1];

# Reactant expects runtime arrays as `RArray`s. We stage the RF samples,
# phantom, and initial controls; the fixed target ``\boldsymbol{b}`` remains constant.

rf_samples_ra = Reactant.to_rarray(p.seq_template.RF[1].A)
obj_ra = Reactant.to_rarray(p.obj)
x0 = Reactant.to_rarray(zeros(ComplexF64, num_nodes));

rf_event = p.seq_template.RF[1] #hide
rf_events = Matrix{RF}(p.seq_template.RF) #hide
rf_events[1] = RF(rf_samples_ra, rf_event.T, rf_event.Δf, rf_event.delay, rf_event.center, rf_event.ϕ, rf_event.use, Val(:preserve)) #hide
seq_ra = Sequence(p.seq_template.GR, rf_events, p.seq_template.ADC, p.seq_template.DUR, p.seq_template.EXT, p.seq_template.DEF) #hide
p_ra = merge(p, (; seq_template=seq_ra, obj=obj_ra)); #hide

# ## Putting everything together
#
# Starting from complex zeros, we form one fixed ``p`` for each of Bloch,
# Magnus2, and Magnus4. A fixed-step gradient descent updates the nodes according to
#
# ```math
# \boldsymbol{x}_{k+1}
# = \boldsymbol{x}_k - \eta\,\nabla_{\boldsymbol{x}}\mathcal{L}_p(\boldsymbol{x}_k).
# ```
#
# We use fifteen steps with ``\eta=1250``. Larger step sizes destabilize the
# lower-order simulation methods. Only the simulation method changes between runs.

# This simple loop can be replaced by any optimizer.
num_iterations = 15
η = 1250.0
function gradient_descent(∇𝓛p, x, p, num_iterations, η)
    for _ in 1:num_iterations
        x = x .- η .* ∇𝓛p(x, p)
    end
    return x
end;

# Test the optimization with three simulation methods.
method_names = ("Bloch-Opt", "Magnus2-Opt", "Magnus4-Opt")
sim_methods = (Bloch(), BlochMagnus2(), BlochMagnus4())
original_stderr = stderr #hide
redirect_stderr(devnull) #hide
optimized_controls = map(sim_methods) do sim_method
    sim_params = merge(p_ra.sim_params, Dict{String,Any}("sim_method" => sim_method))
    p = merge(p_ra, (; sim_params))
    compiled_∇𝓛p = Reactant.@allowscalar Reactant.compile(∇𝓛p, (x0, p))
    gradient_descent(compiled_∇𝓛p, x0, p, num_iterations, η)
end;
redirect_stderr(original_stderr); #hide

# The final figure compares the original Hamming-windowed sinc pulse and the
# three optimized pulses against a common ground-truth magnetization. Each pulse is
# re-simulated with Magnus6 in Float64, using
# ``\Delta t_{\mathrm{rf}}=1\,\mu\mathrm{s}`` and four times the spatial
# resolution of the optimization grid. The reported NRMSE evaluates the full
# complex ``M_{xy}``, not either displayed component alone.

reference_num_spins = 4 * (Nspins - 1) + 1 #hide
reference_z = collect(range(-z_extent, z_extent; length=reference_num_spins)) #hide
reference_b = @. im / sqrt(1 + (abs(reference_z) / z_c)^(2 * butterworth_degree)) #hide
reference_sim_params = merge(p.sim_params, Dict{String,Any}( #hide
    "sim_method" => BlochMagnus6(), #hide
    "precision" => "f64", #hide
    "Δt_rf" => 1e-6, #hide
)) #hide
reference_p = merge(p, (; #hide
    obj=Phantom(; x=zeros(reference_num_spins), z=reference_z), #hide
    sim_params=reference_sim_params, #hide
)) #hide
optimized_reference_profiles = map(optimized_controls) do controls #hide
    𝓐p(Array(controls), reference_p) #hide
end #hide
sinc_profile = simulate(reference_p.obj, reference_p.seq_template, reference_p.sys; sim_params=reference_p.sim_params, verbose=false).xy #hide
comparison_names = ("Sinc", method_names...) #hide
reference_profiles = (sinc_profile, optimized_reference_profiles...) #hide
ground_truth_nrmse_percent = NamedTuple{Symbol.(comparison_names)}(map(reference_profiles) do profile #hide
    100 * sqrt(sum(abs2, profile .- reference_b) / sum(abs2, reference_b)) #hide
end) #hide
display_times = range(first(node_times), last(node_times); length=10 * num_nodes) #hide
optimized_pulses = map(optimized_controls) do controls #hide
    cubic_interpolate_samples((t=node_times, A=Array(controls)), display_times) #hide
end #hide
sinc_display_indices = round.(Int, range(1, length(rf.A); length=length(display_times))) #hide
comparison_times = (sample_times[sinc_display_indices], ntuple(_ -> display_times, length(method_names))...) #hide
comparison_pulses = (1e6 .* rf.A[sinc_display_indices], optimized_pulses...) #hide
target_Mx_trace = scatter(x=reference_z .* 1e3, y=zeros(reference_num_spins), name="Target", legendgroup="Target", legendrank=0, line=attr(color="#AB63FA", dash="dash"), xaxis="x2", yaxis="y2") #hide
target_My_trace = scatter(x=reference_z .* 1e3, y=imag.(reference_b), name="Target", legendgroup="Target", line=attr(color="#AB63FA", dash="dash"), showlegend=false, xaxis="x3", yaxis="y3") #hide
pulse_colors = ("#FFA15A", "#636EFA", "#00CC96", "#EF553B") #hide
node_colors = ("#3B4CC0", "#009E73", "#B33A2B") #hide
pulse_traces = map(eachindex(comparison_names)) do i #hide
    scatter(x=comparison_times[i] .* 1e3, y=real.(comparison_pulses[i]), name=comparison_names[i], legendgroup=comparison_names[i], legendrank=i, line=attr(color=pulse_colors[i]), xaxis="x", yaxis="y") #hide
end #hide
node_traces = map(eachindex(method_names)) do i #hide
    scatter(x=node_times .* 1e3, y=real.(Array(optimized_controls[i])), name=method_names[i], legendgroup=method_names[i], mode="markers", marker=attr(color=node_colors[i], size=6, symbol="circle"), showlegend=false, xaxis="x", yaxis="y") #hide
end #hide
Mx_traces = map(eachindex(comparison_names)) do i #hide
    scatter(x=reference_z .* 1e3, y=real.(reference_profiles[i]), name=comparison_names[i], legendgroup=comparison_names[i], line=attr(color=pulse_colors[i]), showlegend=false, xaxis="x2", yaxis="y2") #hide
end #hide
My_traces = map(eachindex(comparison_names)) do i #hide
    scatter(x=reference_z .* 1e3, y=imag.(reference_profiles[i]), name=comparison_names[i], legendgroup=comparison_names[i], line=attr(color=pulse_colors[i]), showlegend=false, xaxis="x3", yaxis="y3") #hide
end #hide
nrmse_labels = join(("<b>NRMSE</b>", comparison_names...), "<br>") #hide
nrmse_percentages = join(vcat("&nbsp;", [@sprintf("%.2f %%", ground_truth_nrmse_percent[i]) for i in eachindex(comparison_names)]), "<br>") #hide
magnetization_range = [-0.3, 1.1] #hide
profile_position_range = [-2.0, 2.0] #hide
profile_plot = Plot( #hide
    vcat(pulse_traces, node_traces, [target_Mx_trace], Mx_traces, [target_My_trace], My_traces), #hide
    Layout( #hide
        xaxis=attr(title="time [ms]", domain=[0.0, 1.0], anchor="y"), #hide
        yaxis=attr(title="Re(B1) [μT]", domain=[0.58, 1.0], anchor="x"), #hide
        xaxis2=attr(title="z [mm]", domain=[0.0, 0.47], anchor="y2", range=profile_position_range), #hide
        yaxis2=attr(domain=[0.0, 0.42], anchor="x2", range=magnetization_range), #hide
        xaxis3=attr(title="z [mm]", domain=[0.53, 1.0], anchor="y3", range=profile_position_range), #hide
        yaxis3=attr(domain=[0.0, 0.42], anchor="x3", range=magnetization_range), #hide
        legend=attr(groupclick="togglegroup"), #hide
        shapes=[attr(type="rect", x0=1.01, x1=1.32, y0=0.26, y1=0.43, xref="paper", yref="paper", line=attr(color="rgba(0,0,0,0.25)", width=1), fillcolor="rgba(255,255,255,0.85)")], #hide
        annotations=[ #hide
            attr(text="<b>RF pulses</b>", x=0.5, y=1.04, xref="x domain", yref="y domain", xanchor="center", yanchor="bottom", showarrow=false), #hide
            attr(text="<b>M<sub>x</sub></b>", x=0.5, y=1.04, xref="x2 domain", yref="y2 domain", xanchor="center", yanchor="bottom", showarrow=false), #hide
            attr(text="<b>M<sub>y</sub></b>", x=0.5, y=1.04, xref="x3 domain", yref="y3 domain", xanchor="center", yanchor="bottom", showarrow=false), #hide
            attr(text=nrmse_labels, x=1.02, y=0.42, xref="paper", yref="paper", xanchor="left", yanchor="top", align="left", showarrow=false), #hide
            attr(text=nrmse_percentages, x=1.31, y=0.42, xref="paper", yref="paper", xanchor="right", yanchor="top", align="right", showarrow=false), #hide
        ], #hide
        height=700, #hide
        margin=attr(r=230, t=60), #hide
    ), #hide
) #hide
#jl display(profile_plot)
#md profile_plot #hide

# The Magnus methods produce more accurate pulses than Bloch's hard-pulse
# approximation because of their better numerical properties. We intentionally
# used a coarse ``\Delta t_{\mathrm{rf}}`` to make the differences more noticeable,
# but this tutorial demonstrates how KomaMRI's differentiable simulations can
# optimize an RF pulse for a desired target profile.

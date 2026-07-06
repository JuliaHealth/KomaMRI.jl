# # Magnus methods

# ```@meta
# CurrentModule = KomaMRICore
# ```

# Magnus methods are RF-excitation solvers that keep each time step as a
# rotation. They are useful when RF, gradients, or off-resonance change during a
# time step: you can either use a larger ``\Delta t`` for speed, or keep the same
# ``\Delta t`` and reduce time-discretization error.

# ```@raw html
# <p align="center">
#   <img width="100%" src="../assets/magnus-methods.svg">
# </p>
# ```

# KomaMRI provides four recommended Magnus variants:

# | Method | Field model inside one time step | Practical use |
# |:---|:---|:---|
# | [`BlochMagnus1()`](@ref BlochMagnusConst1) | piecewise constant | hard-pulse approximation / baseline |
# | [`BlochMagnus2()`](@ref BlochMagnusMid2) | midpoint field | low-cost second-order default |
# | [`BlochMagnus4()`](@ref BlochMagnusGL4) | Gauss-Legendre samples plus commutator correction | high accuracy, especially with `Float32` GPU simulations |
# | [`BlochMagnus6()`](@ref BlochMagnusBGL6) | Gauss-Legendre samples plus higher-order commutators | highest smooth-order method, most useful with `Float64` |

# Use them through `sim_params["sim_method"]`:

# ```julia
# sim_params = KomaMRICore.default_sim_params()
# sim_params["sim_method"] = BlochMagnus2()
# sim_params["Δt_rf"] = 8e-6
# raw = simulate(obj, seq, sys; sim_params)
# ```

# ## Accuracy And Step Size

# The practical benefit is that higher-order Magnus methods can either improve
# RF excitation accuracy at the same ``\Delta t``, or reach similar accuracy with
# a larger ``\Delta t``.
using KomaMRI, PlotlyJS #hide
const HS_DURATION_DOC = 18.3e-3 #hide
const DT_RF_DOC = [300e-6, 250e-6, 200e-6, 150e-6, 100e-6, 50e-6, 20e-6, 10e-6, 5e-6, 2e-6, 1e-6] #hide
magnus_methods_doc() = ( #hide
    ("Magnus1", BlochMagnus1()), #hide
    ("Magnus2", BlochMagnus2()), #hide
    ("Magnus4", BlochMagnus4()), #hide
    ("Magnus6", BlochMagnus6()), #hide
) #hide
function hs_sequence_doc(n=400_001) #hide
    β̂, μ = 4, 6 #hide
    β = 2β̂ / HS_DURATION_DOC #hide
    t = range(-HS_DURATION_DOC / 2, HS_DURATION_DOC / 2; length=n) #hide
    seq = Sequence() #hide
    @addblock seq += RF(13.5e-6 .* sech.(β .* t), HS_DURATION_DOC, -μ * β .* tanh.(β .* t) ./ (2π), 0) #hide
    return seq #hide
end #hide
function phantom_grid_doc(nspins=100) #hide
    return Phantom(; #hide
        x=collect(range(-6e-2, 6e-2; length=nspins)), #hide
        y=zeros(nspins), #hide
        z=zeros(nspins), #hide
        ρ=ones(nspins), #hide
        T1=fill(Inf, nspins), #hide
        T2=fill(Inf, nspins), #hide
        Δw=2π .* collect(range(-2e3, 2e3; length=nspins)), #hide
    ) #hide
end #hide
function sim_params_doc(method, Δt_rf; precision) #hide
    return Dict{String,Any}( #hide
        "gpu" => false, #hide
        "Nthreads" => Threads.nthreads(), #hide
        "return_type" => "state", #hide
        "precision" => precision, #hide
        "sim_method" => method, #hide
        "Δt" => Inf, #hide
        "Δt_rf" => Δt_rf, #hide
    ) #hide
end #hide
nrmse_doc(M, Mref) = sqrt(sum(abs2.(M.xy .- Mref.xy) .+ abs2.(M.z .- Mref.z)) / sum(abs2.(Mref.xy) .+ abs2.(Mref.z))) #hide
function convergence_result_doc() #hide
    seq, obj, sys = hs_sequence_doc(), phantom_grid_doc(), Scanner() #hide
    reference = simulate(obj, seq, sys; sim_params=sim_params_doc(BlochMagnus6(), 1e-6; precision="f64"), verbose=false) #hide
    methods = magnus_methods_doc() #hide
    err64 = [ #hide
        nrmse_doc(simulate(obj, seq, sys; sim_params=sim_params_doc(method, Δt_rf; precision="f64"), verbose=false), reference) #hide
        for (_, method) in methods, Δt_rf in DT_RF_DOC #hide
    ] #hide
    err32 = [ #hide
        nrmse_doc(simulate(obj, seq, sys; sim_params=sim_params_doc(method, Δt_rf; precision="f32"), verbose=false), reference) #hide
        for (_, method) in methods, Δt_rf in DT_RF_DOC #hide
    ] #hide
    return Dict( #hide
        "method_labels" => collect(first.(methods)), #hide
        "dt_rf_us" => 1e6 .* DT_RF_DOC, #hide
        "err64" => err64, #hide
        "err32" => err32, #hide
    ) #hide
end #hide
result = convergence_result_doc() #hide
labels = result["method_labels"] #hide
dt_rf_us = result["dt_rf_us"] #hide
colors = ["#D55E00", "#009E73", "#0072B2", "#CC79A7"] #hide
traces = GenericTrace[] #hide
for (i, (label, color)) in enumerate(zip(labels, colors)) #hide
    f64_dt = label == "Magnus6" ? dt_rf_us[begin:(end - 1)] : dt_rf_us #hide
    f64_err = label == "Magnus6" ? result["err64"][i, begin:(end - 1)] : result["err64"][i, :] #hide
    push!(traces, scatter(; x=f64_dt, y=f64_err, mode="lines+markers", name="$label Float64", line=attr(color=color), marker=attr(color=color))) #hide
    push!(traces, scatter(; x=dt_rf_us, y=result["err32"][i, :], mode="lines+markers", name="$label Float32", line=attr(color=color, dash="dash"), marker=attr(color=color))) #hide
end #hide
p = plot( #hide
    traces, #hide
    Layout(; #hide
        title="Method convergence", #hide
        xaxis=attr(title="Δt (μs)", type="log", autorange="reversed", tickvals=[300, 100, 50, 20, 10, 5, 2, 1]), #hide
        yaxis=attr(title="NRMSE", type="log", range=[-13, -0.4]), #hide
        width=620, #hide
        height=475, #hide
        margin=attr(l=60, r=150, t=45, b=55), #hide
        legend=attr(x=1.02, y=1, xanchor="left", yanchor="top"), #hide
    ), #hide
) #hide
HTML(replace(sprint(show, MIME("text/html"), p), "style=\"display:block;" => "style=\"display:block;margin-left:auto;margin-right:auto;"; count=1)) #hide

# Solid lines are `Float64`; dashed lines are `Float32`. The convergence
# behavior follows the expected order. In `Float64`, the higher-order methods
# keep improving as ``\Delta t`` becomes smaller. In `Float32`, the error
# eventually reaches roundoff limits, so decreasing ``\Delta t`` indefinitely is
# not always useful.

# ## Sampling Inside A Simulation Step

# A Magnus method may request additional waveform evaluations inside a simulation
# interval; those waveform values are combined into the rotation applied for
# that step.
seq = hs_sequence_doc(4001) #hide
methods = magnus_methods_doc() #hide
plots = map(methods) do (_, method) #hide
    params = Dict{String,Any}("sim_method" => method, "Δt" => Inf, "Δt_rf" => 300e-6) #hide
    rule = KomaMRICore.simulation_sampling_rule(method, params) #hide
    plot_seqd(seq; sampling_rule=rule, show_rf_frame=false) #hide
end #hide
p = KomaMRIPlots.PlotlyJS.make_subplots( #hide
    rows=1, #hide
    cols=4, #hide
    subplot_titles=reshape(collect(first.(methods)), 1, :), #hide
    shared_xaxes=true, #hide
    shared_yaxes=true, #hide
    horizontal_spacing=0.04, #hide
) #hide
for (col, panel) in enumerate(plots), trace in panel.plot.data #hide
    KomaMRIPlots.PlotlyJS.add_trace!(p, trace, row=1, col=col) #hide
end #hide
for axis in (:xaxis, :xaxis2, :xaxis3, :xaxis4) #hide
    p.plot.layout[axis][:range] = [7.5, 9.5] #hide
end #hide
for axis in (:yaxis, :yaxis2, :yaxis3, :yaxis4) #hide
    p.plot.layout[axis][:range] = [10, 14] #hide
end #hide
for axis in (:xaxis2, :xaxis3, :xaxis4) #hide
    p.plot.layout[axis][:matches] = "x" #hide
end #hide
for axis in (:yaxis2, :yaxis3, :yaxis4) #hide
    p.plot.layout[axis][:matches] = "y" #hide
end #hide
KomaMRIPlots.PlotlyJS.relayout!( #hide
    p; #hide
    height=360, #hide
    width="100%", #hide
    showlegend=false, #hide
    xaxis_title="time (ms)", #hide
    xaxis2_title="time (ms)", #hide
    xaxis3_title="time (ms)", #hide
    xaxis4_title="time (ms)", #hide
) #hide
p #hide

# Circles mark simulation-step boundaries. Vertical markers mark internal
# evaluation nodes. The node positions are normalized to one simulation interval
# from 0 to 1.
integration_nodes = KomaMRICore.integration_nodes #hide
eval_intervals_per_step = KomaMRICore.eval_intervals_per_step #hide
println("method   integration_nodes (0 to 1)                         intervals") #hide
for (name, method) in methods #hide
    println(rpad(name, 9), rpad(string(integration_nodes(method)), 56), eval_intervals_per_step(method)) #hide
end #hide

# To add a new Magnus method with internal RF samples, define `integration_nodes`,
# set `eval_intervals_per_step` when the default interval count is not sufficient,
# and implement the corresponding rotation vector.

# ## Effective Field Vector

# Let ``\boldsymbol{B}(t)`` be the effective field in the RF rotating frame,
# including RF frequency modulation. In angular-frequency units,

# ```math
# \boldsymbol{\omega}(t) = -\gamma \boldsymbol{B}(t) \:.
# ```

# ```@raw html
# <p align="center">
#   <img width="58%" src="../assets/magnus-effective-field.png">
# </p>
# ```

# Here ``\boldsymbol{B}(t)`` combines the transverse RF components
# ``B_{1,x}`` and ``B_{1,y}`` with the longitudinal frequency-offset term
# ``(\omega - \omega_{\mathrm{RF}})/\gamma``.

# For one excitation step, the Magnus expansion approximates the time-ordered
# field integral by a skew-symmetric generator:

# ```math
# \boldsymbol{M}_{n+1} =
# \exp\left([\boldsymbol{\theta}]_\times\right) \boldsymbol{M}_n \:,
# ```

# where ``[\cdot]_\times`` is the cross-product matrix. The implementation
# computes ``\boldsymbol{\theta}`` directly. The Magnus rotation vector produces
# a rotation around ``\hat{\boldsymbol{n}}`` with angle ``\theta``:

# ```math
# \exp\left([\boldsymbol{\theta}]_\times\right)
# =
# R_{\hat{\boldsymbol{n}}}(\theta),
# \quad
# \hat{\boldsymbol{n}} =
# \frac{\boldsymbol{\theta}}{\lVert\boldsymbol{\theta}\rVert},
# \quad
# \theta = \lVert\boldsymbol{\theta}\rVert \:.
# ```

# ## Implemented Methods

# The constant, linear, and quadratic methods model
# ``\boldsymbol{\omega}(\tau)`` inside the step.

# | Method | Waveform model | Magnus terms | Smooth order | Additional in-interval evaluations |
# |:---|:---|:---|---:|---:|
# | [`BlochMagnusConst1`](@ref) | constant | ``\theta_1`` | ``\mathcal{O}(\Delta t)`` | 0 |
# | [`BlochMagnusLin2`](@ref) | linear | ``\theta_1`` | ``\mathcal{O}(\Delta t^2)`` | 0 |
# | [`BlochMagnusMid2`](@ref) | midpoint | ``\theta_1`` | ``\mathcal{O}(\Delta t^2)`` | 1 |
# | [`BlochMagnusLinComm2`](@ref) | linear | ``\theta_1+\theta_2`` | ``\mathcal{O}(\Delta t^2)`` | 0 |
# | [`BlochMagnusQuad2`](@ref) | quadratic | ``\theta_1`` | ``\mathcal{O}(\Delta t^2)`` | 1 |
# | [`BlochMagnusQuad4`](@ref) | quadratic | ``\theta_1+\theta_2`` | ``\mathcal{O}(\Delta t^4)`` | 1 |

# The Gauss-Legendre methods approximate the Magnus integrals directly by
# evaluating ``\boldsymbol{\omega}(\tau)`` at quadrature nodes.

# | Method | Integral approximation | Magnus terms | Smooth order | Additional in-interval evaluations |
# |:---|:---|:---|---:|---:|
# | [`BlochMagnusGL2`](@ref) | Gauss-Legendre | ``\theta_1`` | ``\mathcal{O}(\Delta t^2)`` | 2 |
# | [`BlochMagnusGL4`](@ref) | Gauss-Legendre | ``\theta_1+\theta_2`` | ``\mathcal{O}(\Delta t^4)`` | 2 |
# | [`BlochMagnusBGL4`](@ref) | Blanes Gauss-Legendre | ``\theta_1+\theta_2`` | ``\mathcal{O}(\Delta t^4)`` | 3 |
# | [`BlochMagnusBGL6`](@ref) | Blanes Gauss-Legendre | ``\theta_1+\theta_2+\theta_{34}`` | ``\mathcal{O}(\Delta t^6)`` | 3 |

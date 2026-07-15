using KomaMRI, PlotlyBase #hide
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
p = Plot( #hide
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

seq = hs_sequence_doc(4001) #hide
methods = magnus_methods_doc() #hide
plots = map(methods) do (_, method) #hide
    params = Dict{String,Any}("sim_method" => method, "Δt" => Inf, "Δt_rf" => 300e-6) #hide
    rule = KomaMRICore.simulation_sampling_rule(method, params) #hide
    plot_seqd(seq; sampling_rule=rule, show_rf_frame=false) #hide
end #hide
p = Plot(Layout(Subplots(; #hide
    rows=1, #hide
    cols=4, #hide
    subplot_titles=reshape(collect(first.(methods)), 1, :), #hide
    shared_xaxes=true, #hide
    shared_yaxes=true, #hide
    horizontal_spacing=0.04, #hide
))) #hide
for (col, panel) in enumerate(plots), trace in panel.data #hide
    add_trace!(p, trace; row=1, col=col) #hide
end #hide
for axis in (:xaxis, :xaxis2, :xaxis3, :xaxis4) #hide
    p.layout[axis][:range] = [7.5, 9.5] #hide
end #hide
for axis in (:yaxis, :yaxis2, :yaxis3, :yaxis4) #hide
    p.layout[axis][:range] = [10, 14] #hide
end #hide
for axis in (:xaxis2, :xaxis3, :xaxis4) #hide
    p.layout[axis][:matches] = "x" #hide
end #hide
for axis in (:yaxis2, :yaxis3, :yaxis4) #hide
    p.layout[axis][:matches] = "y" #hide
end #hide
relayout!( #hide
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

integration_nodes = KomaMRICore.integration_nodes #hide
eval_intervals_per_step = KomaMRICore.eval_intervals_per_step #hide
println("method   integration_nodes (0 to 1)                         intervals") #hide
for (name, method) in methods #hide
    println(rpad(name, 9), rpad(string(integration_nodes(method)), 56), eval_intervals_per_step(method)) #hide
end #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

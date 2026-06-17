# # Adiabatic RF Pulse

using KomaMRI #hide
sys = Scanner(); #hide

# In this tutorial, we will build a hyperbolic-secant (HS) adiabatic inversion
# pulse. The key point is that KomaMRI can keep the RF frequency modulation
# explicit, then simulate the RF-frame dynamics and the resulting
# ``B_0``/``B_1`` robustness.

# ## Defining a frequency-modulated pulse
#
# > Parameters follow the default ``\mathrm{HS}_{4,6}`` pulse in Figure 1(c) of this
# > [JMRI adiabatic inversion example](https://doi.org/10.1002/jmri.26021).

b1max = 13.5e-6
Trf = 18.3e-3
β̂ = 4
μ = 6
β = 2 * β̂ / Trf;

# For an HS pulse, the hardest point for adiabatic following is near the
# center of the sweep, where ``B_1`` is maximal and
# ``|\mathrm{d}\Delta\omega/\mathrm{d}t|=\mu\beta^2``. Requiring the RF precession rate to
# dominate that sweep rate gives the threshold used in the paper:
#
# ```math
# \frac{(\gamma_\mathrm{rad} B_1)^2}{\mu\beta^2} \ge 1
# \quad\Rightarrow\quad
# B_1 \geq \frac{\sqrt{\mu}\beta}{\gamma_\mathrm{rad}}.
# ```

b1_threshold = β * sqrt(μ) / (2π * γ)
b1max > b1_threshold

# First, we define the RF amplitude as a hyperbolic secant and the frequency
# sweep as a hyperbolic tangent.

t = range(-Trf / 2, Trf / 2, 201)
B1 = b1max .* sech.(β .* t)
Δf = -μ * β .* tanh.(β .* t) ./ (2π);

# The ``\Delta f`` argument of `RF` can be a scalar, for a constant RF offset
# such as a slice offset, or a waveform, as in this frequency-modulated pulse.

f = range(-2e3, 2e3, 161) |> collect
seq = Sequence()
@addblock seq += RF(B1, Trf, Δf, 0);

# This representation is native to KomaMRI. Pulseq RF events do not store
# frequency-modulation waveforms, so exporting this sequence with `write_seq`
# requires first converting the sweep into RF phase samples. KomaMRI may do this
# conversion automatically in the future.

# ## Plotting an adiabatic pulse
#
# `plot_seq` can show both views. The default plot keeps ``\Delta f(t)``
# explicit; the `freq_in_phase` keyword shows the same pulse after moving the
# frequency sweep into the RF phase.

using PlotlyJS #hide
function show_only_traces!(p, names) #hide
    for trace in p.plot.data #hide
        name = get(trace.fields, :name, "") #hide
        trace.fields[:showlegend] = false #hide
        trace.fields[:visible] = name in names #hide
    end #hide
    return p #hide
end #hide
function scale_trace_y!(p, name, scale) #hide
    for trace in p.plot.data #hide
        get(trace.fields, :name, "") == name || continue #hide
        y = trace.fields[:y] #hide
        trace.fields[:customdata] = y #hide
        trace.fields[:y] = scale .* y #hide
        trace.fields[:hovertemplate] = "(%{x:.4f} ms, Δf_FM: %{customdata:.4f} kHz)" #hide
    end #hide
    return p #hide
end #hide
p_freq = plot_seq(seq; max_rf_samples=Inf, slider=false, height=360, title="Frequency-modulated RF", showlegend=false)
p_phase = plot_seq(seq; freq_in_phase=true, max_rf_samples=Inf, slider=false, height=360, title="Phase-modulated RF", showlegend=false)
show_only_traces!(p_freq, ("|B1|_AM", "Δf_FM")) #hide
show_only_traces!(p_phase, ("|B1|_AM", "∠B1_AM")) #hide
scale_trace_y!(p_freq, "Δf_FM", 9) #hide
p_rf = [p_freq p_phase] #hide
relayout!(p_rf, height=380, margin=attr(t=56)) #hide
p_rf #hide
#jl display(p_rf)

# Keeping the frequency sweep explicit is useful because KomaMRI simulates a
# frequency-modulated RF pulse in the RF rotating frame. If the same pulse is
# first converted to sampled phase modulation, the simulation is more sensitive
# to the RF time sampling.

# ## Comparing RF and rotating frames
#
# Using a callback, we can record the magnetization during the RF pulse.

trajectory = NamedTuple[]
call_every_N_blocks = 1

record_traj = Callback(
    call_every_N_blocks,
    (progress_info, sim_blocks_info, device_data, sim_params) -> begin
        j = last(sim_blocks_info.parts[progress_info.block])
        push!(trajectory, (;
            Mxy=device_data.Xt.xy[1], Mz=device_data.Xt.z[1],
            ψ=device_data.seqd.ψ[j], B1=device_data.seqd.B1[j], Δf=device_data.seqd.Δf[j],
        ))
    end,
)

sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "state"
sim_params["max_rf_block_length"] = 1; # very inefficient; just for plots
obj0 = Phantom(; x=[0.0], Δw=[0.0]);
simulate(obj0, seq, sys; sim_params, callbacks=(record_traj,), verbose=false);

function adiabatic_frame(p) #hide
    ωeff = (-real(p.B1), -imag(p.B1), p.Δf / γ) #hide
    ω̂rf = ωeff ./ sqrt(sum(abs2, ωeff)) #hide
    Mxy_rf = p.Mxy * cis(-p.ψ) #hide
    ω̂xy_rot = complex(ω̂rf[1], ω̂rf[2]) * cis(p.ψ) #hide
    return (; #hide
        Mrf=(real(Mxy_rf), imag(Mxy_rf), p.Mz), #hide
        Mrot=(real(p.Mxy), imag(p.Mxy), p.Mz), #hide
        ω̂rf, #hide
        ω̂rot=(real(ω̂xy_rot), imag(ω̂xy_rot), ω̂rf[3]), #hide
    ) #hide
end #hide
trajectory_frames = adiabatic_frame.(trajectory) #hide
xyz(points) = ntuple(i -> getindex.(points, i), 3) #hide
Mrf_xyz = xyz(getproperty.(trajectory_frames, :Mrf)) #hide
Mrot_xyz = xyz(getproperty.(trajectory_frames, :Mrot)) #hide
ω̂rf_xyz = xyz(getproperty.(trajectory_frames, :ω̂rf)) #hide
ω̂rot_xyz = xyz(getproperty.(trajectory_frames, :ω̂rot)) #hide
anim_idx = unique(round.(Int, range(1, length(trajectory_frames), 36))) #hide
function sphere_mesh(; nθ=36, nφ=18) #hide
    θ = range(0, 2π, nθ) #hide
    φ = range(0, π, nφ) #hide
    return ( #hide
        x=[sin(p) * cos(t) for p in φ, t in θ], #hide
        y=[sin(p) * sin(t) for p in φ, t in θ], #hide
        z=[cos(p) for p in φ, t in θ], #hide
    ) #hide
end #hide
sphere = sphere_mesh() #hide
function bloch_traces(i, M, ω̂; scene="scene", showlegend=true) #hide
    path_trace = scatter3d(; #hide
        x=M[1][1:i], y=M[2][1:i], z=M[3][1:i], #hide
        mode="lines", name="M path", scene=scene, showlegend=showlegend, #hide
        line=attr(color="#111827", width=6), #hide
    ) #hide
    magnetization = scatter3d(; #hide
        x=[0, M[1][i]], y=[0, M[2][i]], z=[0, M[3][i]], #hide
        mode="lines+markers", name="M", scene=scene, showlegend=showlegend, #hide
        line=attr(color="#dc2626", width=7), #hide
        marker=attr(color="#dc2626", size=4), #hide
    ) #hide
    effective_field = scatter3d(; #hide
        x=[0, ω̂[1][i]], y=[0, ω̂[2][i]], z=[0, ω̂[3][i]], #hide
        mode="lines+markers", name="ωeff", scene=scene, showlegend=showlegend, #hide
        line=attr(color="#2563eb", width=7), #hide
        marker=attr(color="#2563eb", size=4), #hide
    ) #hide
    return [ #hide
        path_trace, #hide
        magnetization, #hide
        effective_field, #hide
    ] #hide
end #hide
bloch_sphere(scene) = [surface(; x=sphere.x, y=sphere.y, z=sphere.z, opacity=0.14, showscale=false, name="Bloch sphere", scene=scene, colorscale=[[0, "#dbe4ef"], [1, "#dbe4ef"]])] #hide
function animation_frame(i) #hide
    traces = [ #hide
        bloch_traces(i, Mrf_xyz, ω̂rf_xyz)..., #hide
        bloch_traces(i, Mrot_xyz, ω̂rot_xyz; scene="scene2", showlegend=false)..., #hide
    ] #hide
    return frame(; name=string(i), data=traces, traces=[1, 2, 3, 5, 6, 7]) #hide
end #hide
function bloch_scene(domain) #hide
    axis(range) = attr(; #hide
        range=range, title=attr(text=""), ticks="", #hide
        showgrid=false, showbackground=false, showticklabels=false, zeroline=false, #hide
    ) #hide
    return attr(; #hide
        domain=domain, aspectmode="cube", #hide
        xaxis=axis([-1, 1]), yaxis=axis([-1, 1]), zaxis=axis([-1, 1]), #hide
        camera=attr(eye=attr(x=1.31, y=1.31, z=0.75), up=attr(x=0, y=0, z=1)), #hide
    ) #hide
end #hide
base_traces = [ #hide
    bloch_sphere("scene")..., #hide
    bloch_traces(first(anim_idx), Mrf_xyz, ω̂rf_xyz)..., #hide
    bloch_sphere("scene2")..., #hide
    bloch_traces(first(anim_idx), Mrot_xyz, ω̂rot_xyz; scene="scene2", showlegend=false)..., #hide
] #hide
frames = animation_frame.(anim_idx) #hide
p_bloch = plot( #hide
    base_traces, #hide
    Layout(; #hide
        height=340, margin=attr(l=0, r=0, t=32, b=0), #hide
        scene=bloch_scene(attr(x=[0.0, 0.48], y=[0.0, 0.93])), #hide
        scene2=bloch_scene(attr(x=[0.52, 1.0], y=[0.0, 0.93])), #hide
        annotations=[ #hide
            attr(text="RF frame", x=0.24, y=0.98, xref="paper", yref="paper", xanchor="center", showarrow=false, font=attr(size=18)), #hide
            attr(text="Rotating frame", x=0.76, y=0.98, xref="paper", yref="paper", xanchor="center", showarrow=false, font=attr(size=18)), #hide
        ], #hide
        updatemenus=[ #hide
            attr(; #hide
                type="buttons", direction="right", x=0, y=1.0, #hide
                xanchor="left", yanchor="top", #hide
                buttons=[attr(label="Play", method="animate", args=[nothing, attr(frame=attr(duration=80, redraw=true), fromcurrent=true)])], #hide
            ), #hide
        ], #hide
    ), #hide
    frames, #hide
) #hide
#jl display(p_bloch)

# The blue vector is the normalized rotation axis
# ``\hat{\boldsymbol{\omega}}_\mathrm{eff}``, where
# ``\boldsymbol{\omega}_\mathrm{eff} = -\gamma \mathbf{B}_\mathrm{eff}``. With this sign convention,
# the magnetization precesses right-handed around the plotted axis.

# ## HS adiabatic pulse B0 and B1 robustness
#
# To showcase the off-resonance and ``B_1`` robustness of this type of pulse, we
# can show its effect in a heatmap with ``B_0 \in [-2, 2]\,\mathrm{kHz}`` and
# ``B_1 \in [0, 16]\,\mu\mathrm{T}``:

obj = Phantom(; x=zeros(length(f)), Δw=2π .* f);
b1_scales = range(0.05, 1.2, 47) |> collect;
sim_params = Dict{String, Any}("return_type" => "state");

Mz = map(b1_scales) do scale
    seq_scale = Sequence()
    @addblock seq_scale += RF(scale .* B1, Trf, Δf, 0)
    simulate(obj, seq_scale, sys; sim_params, verbose=false).z
end;

b1_axis = 1e6 .* b1max .* b1_scales #hide
threshold = 1e6 * b1_threshold #hide
Mz_map = round.(reduce(hcat, Mz); digits=3) #hide
p_hs = plot( #hide
    [ #hide
        heatmap(; x=f, y=b1_axis, z=permutedims(Mz_map), zmin=-1, zmax=1, colorscale="RdBu", colorbar=attr(title="Mz")), #hide
        scatter(; x=f, y=fill(threshold, length(f)), name="threshold", mode="lines", hoverinfo="skip", line=attr(color="#a62023", dash="dash", width=3)), #hide
    ], #hide
    Layout(; title="HS<sub>4,6</sub> inversion profile", xaxis_title="Off-resonance [Hz]", yaxis_title="B1,max [μT]", height=420, margin=attr(t=64)), #hide
) #hide
p_hs #hide
#jl display(p_hs)

# The dashed line is the analytic adiabatic threshold for this HS pulse, where we
# can see that the inversion is achieved after the threshold is surpassed.

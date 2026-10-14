"""
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)

Define colors for dark or light mode.

# Arguments
- `darkmode`: (`::Bool`) boolean that selects dark or light mode

# Returns
- `bgcolor`: (`::String`) backgound color
- `text_color`: (`::String`) text color
- `plot_bgcolor`: (`::String`) color background for the plots
- `grid_color`: (`::String`) color of the grids
- `sep_color`: (`::String`) color of separator lines
"""
function theme_chooser(darkmode)
    if darkmode
        bgcolor = "rgb(13,16,17)"
        text_color = "gray"
        plot_bgcolor = "rgb(22,26,29)" #"rgb(33,37,41)"
        grid_color = "rgb(40,52,66)" #rgb(40,40,40)
        sep_color = "white"
    else
        bgcolor = "white"
        text_color = "gray"#"rgb(49,70,101)"
        plot_bgcolor = "rgb(229,236,246)"
        grid_color = "white"
        sep_color = "black"
    end
    return bgcolor, text_color, plot_bgcolor, grid_color, sep_color
end

const PHASE_COLORSCALE = [
    (0.0, "#d64d4d"), (1/6, "#c6af42"), (2/6, "#44b079"),
    (3/6, "#469fc2"), (4/6, "#7861c5"), (5/6, "#bd59a2"), (1.0, "#d64d4d"),
]

const KSPACE_COLORSCALE = let
    positions = range(0, 1; length=10)
    hue = linear_interpolation(positions, [floor(Int, (1-t)*255) for t in positions])
    # Match the legacy WebGL renderer, which parsed the HSV strings as HSL.
    [(t, "#" * PlotlyBase.Colors.hex(PlotlyBase.Colors.HSL(hue(t), 1, 0.5)))
        for t in range(0, 1; length=256)]
end

function generate_seq_time_layout_config(
    title, width, height, range, slider, show_seq_blocks, darkmode; T0, label_to_show=(), non_label_count=8
    )

    num_labels = length(label_to_show)
    # For dropdown, only update label traces, leave non-label traces unchanged
    # Plotly uses 0-based indices: label traces start at non_label_count (0-based)
    label_indices = [i-1 for i in (non_label_count+1):(non_label_count + num_labels)]
    buttons = [
        attr(
            label = "Labels (none)",
            method = "restyle",
            args = [attr(visible = fill(false, num_labels)), label_indices]
        )
    ]
    # Add one button per label fieldname (only one label trace visible at a time)
    append!(buttons, [
        attr(
            label = string(sym),
            method = "restyle",
            args = [attr(visible = [j == i for j in 1:num_labels]), label_indices]
        ) for (i, sym) in enumerate(label_to_show)
    ])


    #LAYOUT
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
    l = Layout(;
        title=title,
        hovermode="closest",
        xaxis_title="",
        modebar=attr(;
            orientation="h",
            yanchor="bottom",
            xanchor="right",
            y=1,
            x=0,
            bgcolor=bgcolor,
            color=text_color,
            activecolor=plot_bgcolor,
        ),
        legend=attr(; orientation="h", yref="container", yanchor="top", xanchor="left", y=1, x=0),
        updatemenus=isempty(label_to_show) ? [] : [
            attr(
                type="dropdown",
                yref="paper",
                xref="paper",
                y=0.99,
                x=0.99,
                yanchor="top",
                xanchor="right",
                align="middle",
                orientation="h",
                bgcolor="white",
                color=text_color,
                buttons=buttons,
            )
        ],
        plot_bgcolor=plot_bgcolor,
        paper_bgcolor=bgcolor,
        xaxis_gridcolor=grid_color,
        yaxis_gridcolor=grid_color,
        xaxis_zerolinecolor=grid_color,
        yaxis_zerolinecolor=grid_color,
        font_color=text_color,
        yaxis_fixedrange=false,
        yaxis_automargin=false,
        xaxis=attr(;
            automargin=false,
            ticksuffix=" ms",
            domain=range[:],
            range=range[:],
            rangeslider=attr(; visible=slider),
            rangeselector=attr(;
                buttons=[
                    attr(; count=1, label="1m", step=10, stepmode="backward"),
                    attr(; step="all"),
                ],
            ),
        ),
        margin=attr(; t=6, l=30, r=6, b=24),
    )
    if show_seq_blocks
        l.xaxis["tickvals"] = T0 * 1e3
    end
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
    #CONFIG
    config = PlotConfig(;
        displaylogo=false,
        toImageButtonOptions=attr(;
            format="svg", # one of png, svg, jpeg, webp
        ).fields,
        modeBarButtonsToRemove=[
            "zoom",
            "select2d",
            "lasso2d",
            "autoScale",
            "resetScale2d",
            "pan",
            "tableRotation",
            "resetCameraLastSave",
            "zoomIn",
            "zoomOut",
        ],
    )

    return l, config
end

"""
    c_map_interp = interp_map(c_map, t_interp)

Interpolates a color map. This is used for plotting the kspace (refer to
[`plot_kspace`](@ref)).

# Arguments
- `c_map`: (`::Vector{Vector{Any}}`)color map. Every element of this vector has a
    vector with a number between 0-1 in its first element and a color string in its second
    element. It serves as a reference to create a color map with more elements
- `t_interp`: (`::Vector{Float64}`) the vector with values between 0-1 that are the
    reference for interpolate the color map with more elements

# Returns
- `c_map_interp`: (`::Vector{String}`) vector with color strings with interpolated
    values
"""
function interp_map(c_map, t_interp)
    idx = [c[1] for c in c_map]
    R = [parse.(Int, split(c[2][5:(end - 1)], ", "))[1] for c in c_map]
    G = [parse.(Int, split(c[2][5:(end - 1)], ", "))[2] for c in c_map]
    B = [parse.(Int, split(c[2][5:(end - 1)], ", "))[3] for c in c_map]
    r = linear_interpolation(idx, R)(t_interp)
    g = linear_interpolation(idx, G)(t_interp)
    b = linear_interpolation(idx, B)(t_interp)
    c_map_interp = ["hsv($r, $g, $b)" for (r, g, b) in zip(r, g, b)]
    return c_map_interp
end

_trigger_delays(::Extension) = ()
_trigger_delays(trigger::Trigger) = (trigger.delay,)

function trigger_name(trigger)
    trigger.type == 1 && return "Digital output"
    trigger.type == 2 && return "Input trigger"
    return "Trigger (type $(trigger.type))"
end

function trigger_markers(seq, starts)
    return [(; time=(starts[b] + delay) * 1e3,
        label="$(trigger_name(extension)) · channel $(extension.channel)")
        for b in eachindex(seq.DUR) for extension in seq.EXT[b]
        for delay in _trigger_delays(extension)]
end

_ecg_gaussian(t, center, width) = exp(-0.5 * ((t - center) / width)^2)

function _ecg_r_peaks(signal::CardiacSignal, t_start, t_end)
    if isnothing(signal.period)
        return signal.r_peaks[(t_start - 0.5 .<= signal.r_peaks) .& (signal.r_peaks .<= t_end + 0.5)]
    end

    reference = first(signal.r_peaks)
    first_period = floor(Int, (t_start - 0.5 - reference) / signal.period)
    last_period = ceil(Int, (t_end + 0.5 - reference) / signal.period)
    return reference .+ (first_period:last_period) .* signal.period
end

function _ecg_waveform(t, signal::CardiacSignal)
    ecg = zeros(length(t))
    # Beyond eight Gaussian widths, each component is below exp(-32) of its peak.
    tail_widths = 8
    beat_start = -0.200 - tail_widths * 0.035
    beat_end = 0.280 + tail_widths * 0.070
    for r_peak in _ecg_r_peaks(signal, first(t), last(t))
        first_sample = searchsortedfirst(t, r_peak + beat_start)
        last_sample = searchsortedlast(t, r_peak + beat_end)
        @inbounds for i in first_sample:last_sample
            relative_time = t[i] - r_peak
            ecg[i] +=
                0.12 * _ecg_gaussian(relative_time, -0.200, 0.035) -
                0.15 * _ecg_gaussian(relative_time, -0.040, 0.012) +
                1.00 * _ecg_gaussian(relative_time, 0.000, 0.010) -
                0.25 * _ecg_gaussian(relative_time, 0.035, 0.014) +
                0.30 * _ecg_gaussian(relative_time, 0.280, 0.070)
        end
    end
    return ecg
end

_add_physio!(p, layout, scatter_fun, seq, ::NoPhysioSignal, xaxis) = nothing

function _add_physio!(p, layout, scatter_fun, seq, signal::CardiacSignal, xaxis)
    seq_duration = dur(seq)
    r_peaks = filter(peak -> 0.0 <= peak <= seq_duration, _ecg_r_peaks(signal, 0.0, seq_duration))
    interval_edges = unique([0.0; r_peaks; seq_duration])
    samples_per_rr = 400
    ecg_time = [
        interval_start + local_time
        for (interval_start, interval_end) in zip(interval_edges, @view interval_edges[2:end])
        for local_time in LinRange(0.0, interval_end - interval_start, samples_per_rr)
    ]
    push!(p, scatter_fun(;
        x=ecg_time * 1e3,
        y=_ecg_waveform(ecg_time, signal),
        name="ECG",
        hovertemplate="(%{x:.4f} ms, %{y:.3f})<extra>ECG</extra>",
        xaxis=xaxis,
        yaxis="y2",
        showlegend=false,
        line=attr(; color="red"),
    ))
    layout.yaxis["domain"] = [0.12, 1.0]
    layout.xaxis["anchor"] = "y2"
    layout.xaxis[:rangeslider][:visible] = false
    layout.yaxis2 = attr(;
        domain=[0.0, 0.09],
        anchor=xaxis,
        fixedrange=true,
        showgrid=false,
        showticklabels=false,
        zeroline=false,
        title=attr(; text="Physio", font=attr(; size=11), standoff=0),
    )
    return nothing
end

function stack_plot_samples(samples)
    T = mapreduce(eltype, promote_type, samples)
    values = Vector{Union{Missing,T}}(undef, sum(length, samples) + length(samples))
    offset = 1
    for sample in samples
        copyto!(values, offset, sample, 1, length(sample))
        offset += length(sample)
        values[offset] = missing
        offset += 1
    end
    return values
end







# Nonadaptive Plotly figure builders, also reused by adaptive time plots.
const RF_CENTER_MARKER = (color="#FF0000", symbol="x")

function plot_seq_nonadaptive(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=false,
    show_seq_blocks=false,
    darkmode=false,
    range=[],
    title="",
    xaxis="x",
    yaxis="y",
    showlegend=true,
    freq_in_phase=false,
    show_rf_frame=false,
    # Performance related
    gl=false,
    max_rf_samples=100,
    show_adc=false,
    physio=NoPhysioSignal(),
)

    seq = resolve_triggers(seq, physio)
    # Aux functions
    scatter_fun = gl ? scattergl : scatter
    usrf(x) = length(x) > max_rf_samples ? ([@view x[1]; @view x[2:(length(x)÷max_rf_samples):end-1]; @view x[end]]) : x
    usadc(x; ampl_edge=1.0) = show_adc || isempty(x) ? x : [ampl_edge * first(x); 1.0 * first(x); 1.0 * last(x); ampl_edge * last(x)]
    # Sample blocks locally so their start times are not recomputed for the full sequence.
    T0 = get_block_start_times(seq)
    seq_samples = [get_samples(block; freq_in_phase) for block in seq]
    stack_samples(name; amp=identity, time=identity) = (
        A=stack_plot_samples([amp(getproperty(block, name).A) for block in seq_samples]),
        t=stack_plot_samples([time(T0[i] .+ getproperty(block, name).t)
            for (i, block) in enumerate(seq_samples)]),
    )
    active(values, enabled) = enabled ? values : fill(missing, length(values))
    triggers = trigger_markers(seq, T0)
    trigger_times = getproperty.(triggers, :time)
    # Get center times
    center_times = similar(T0, 0)
    center_values = ComplexF64[]
    for (i, b) in enumerate(seq)
        if is_RF_on(b)
            center_time = b.RF[1].delay + b.RF[1].center
            push!(center_times, T0[i] + center_time)
            push!(center_values, KomaMRIBase.get_rfs(b, [center_time])[1][1])
        end
    end
    gx = stack_samples(:gx)
    gy = stack_samples(:gy)
    gz = stack_samples(:gz)
    rf = (;
        stack_samples(:rf; amp=usrf, time=usrf)...,
        ct=center_times,
        cA=abs.(center_values),
        cϕ=angle.(center_values)
    )
    Δf = stack_samples(:Δf; amp=usrf, time=usrf)
    ψ = show_rf_frame ? stack_samples(:ψ; amp=usrf, time=usrf) : nothing
    adc = stack_samples(:adc; amp=(x -> usadc(x; ampl_edge=0.0)), time=usadc)

    label = get_labels(seq)
    isadc = is_ADC_on.(seq)
    label_symbols = [
        sym for sym in fieldnames(AdcLabels)
        if any(j -> isadc[j] && !iszero(getfield(label[j], sym)), eachindex(label))
    ]

    # Define general params and the vector of plots
    idx = ["Gx" "Gy" "Gz"]
    O = size(seq.RF, 1)
    rf_trace_count = 3 + (freq_in_phase ? 0 : 1) + (!freq_in_phase && show_rf_frame ? 1 : 0)
    adc_idx = 3 + rf_trace_count * O + 1
    p = [scatter_fun() for _ in 1:(adc_idx + length(label_symbols))]

    # For GRADs
    p[1] = scatter_fun(;
        x=gx.t * 1e3,
        y=active(gx.A * 1e3, is_Gx_on(seq)),
        name=idx[1],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
        xaxis=xaxis,
        yaxis=yaxis,
        legendgroup="Gx",
        showlegend=showlegend,
        marker=attr(; color="#636EFA", size=8),
    )
    p[2] = scatter_fun(;
        x=gy.t * 1e3,
        y=active(gy.A * 1e3, is_Gy_on(seq)),
        name=idx[2],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
        xaxis=xaxis,
        yaxis=yaxis,
        legendgroup="Gy",
        showlegend=showlegend,
        marker=attr(; color="#EF553B", size=8),
    )
    p[3] = scatter_fun(;
        x=gz.t * 1e3,
        y=active(gz.A * 1e3, is_Gz_on(seq)),
        name=idx[3],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
        xaxis=xaxis,
        yaxis=yaxis,
        legendgroup="Gz",
        showlegend=showlegend,
        marker=attr(; color="#00CC96", size=8),
    )

    # For RFs
    rf_on = is_RF_on(seq)
    for j in 1:O
        idx_rf = 3 + rf_trace_count * (j - 1)
        rf_wave = rf.A[:, j]
        is_real_rf = all(x -> ismissing(x) || isapprox(imag(x), 0; atol=eps()), rf_wave)
        rf_amp = map(x -> ismissing(x) ? missing : (is_real_rf ? real(x) : abs(x)), rf_wave)
        rf_phase = map(x -> ismissing(x) ? missing : (is_real_rf ? zero(real(x)) : angle(x)), rf_wave)
        # Plot RF
        p[idx_rf + 1] = scatter_fun(;
            x=rf.t * 1e3,
            y=active(rf_amp * 1e6, rf_on),
            name="|B1|_AM",
            hovertemplate="(%{x:.4f} ms, %{y:.2f} μT)",
            xaxis=xaxis,
            yaxis=yaxis,
            legendgroup="|B1|_AM",
            showlegend=showlegend,
            marker=attr(; color="#AB63FA"),
        )
        p[idx_rf + 2] = scatter_fun(;
            x=rf.t * 1e3,
            y=active(rf_phase, rf_on),
            text=ones(size(rf.t)),
            name="∠B1_AM",
            hovertemplate="(%{x:.4f} ms, ∠B1: %{y:.4f} rad)",
            visible="legendonly",
            xaxis=xaxis,
            yaxis=yaxis,
            legendgroup="∠B1_AM",
            showlegend=showlegend,
            marker=attr(; color="#FFA15A"),
        )
        center_idx = idx_rf + 3
        if !freq_in_phase
            p[center_idx] = scatter_fun(;
                x=Δf.t * 1e3,
                y=active(Δf.A[:, j] * 1e-3, rf_on),
                text=ones(size(Δf.t)),
                name="Δf_FM",
                hovertemplate="(%{x:.4f} ms, Δf_FM: %{y:.4f} kHz)",
                visible="legendonly",
                xaxis=xaxis,
                yaxis=yaxis,
                legendgroup="Δf_FM",
                showlegend=showlegend,
                marker=attr(; color="#AB63FA"),
                line=attr(; dash="dot"),
            )
            center_idx += 1
            if show_rf_frame
                p[center_idx] = scatter_fun(;
                    x=ψ.t * 1e3,
                    y=active(ψ.A[:, j], rf_on),
                    text=ones(size(ψ.t)),
                    name="ψ_FM",
                    hovertemplate="(%{x:.4f} ms, ψ_FM: %{y:.4f} rad)",
                    visible="legendonly",
                    xaxis=xaxis,
                    yaxis=yaxis,
                    legendgroup="ψ_FM",
                    showlegend=showlegend,
                    marker=attr(; color="#FF6692"),
                    line=attr(; dash="dot"),
                )
                center_idx += 1
            end
        end
        p[center_idx] = scatter_fun(;
            x=rf.ct * 1e3,
            y=active(rf.cA * 1e6, rf_on),
            text=rf.cϕ,
            name="RF_center",
            hovertemplate="RF center: %{x:.4f} ms<br>|B1|: %{y:.2f} μT<br>∠B1: %{text:.2f} rad<extra></extra>",
            visible="legendonly",
            xaxis=xaxis,
            yaxis=yaxis,
            legendgroup="RF_center",
            showlegend=showlegend,
            mode="markers",
            marker=attr(; RF_CENTER_MARKER...),
        )
    end

    # For ADCs
    p[adc_idx] = scatter_fun(;
        x=adc.t * 1e3,
        y=active(adc.A, is_ADC_on(seq)),
        name="ADC",
        hovertemplate="(%{x:.4f} ms, %{y:i})",
        xaxis=xaxis,
        yaxis=yaxis,
        legendgroup="ADC",
        showlegend=showlegend,
        mode=(show_adc ? "markers" : "line"),
        marker=attr(; color="#19D3F3"),
    )

    #############################
    ###### show label
    ############################
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)

    d = [ seq[i].DUR[1] for i in eachindex(seq.DUR)]
    d2 = [0;d]
    dcum = cumsum(d2)
    t_center = dcum[1:end-1] + d/2
    t_center_adc = t_center[isadc]

    for (i, sym) in enumerate(label_symbols)
        lab_vec = [getfield(label[j],sym) for j in eachindex(label)]
        lab_adc = lab_vec[isadc]

        p[adc_idx + i] = scatter_fun(;
            x= t_center_adc * 1e3,
            y= lab_adc,
            name=string(sym),
            hovertemplate="(%{x:.4f} ms, %{y:i})",
            xaxis=xaxis,
            yaxis=yaxis,
            legendgroup=string(sym),
            showlegend=false,
            mode=("markers"),
            marker=attr(; color=sep_color, symbol="x"),
            meta=attr(; theme="contrast"),
            visible=false,
        )
    end

    ###############################
    ###############################
    # Return the plot
    l, config = generate_seq_time_layout_config(
        title,
        width,
        height,
        range,
        slider,
        show_seq_blocks,
        darkmode;
        T0=get_block_start_times(seq),
        label_to_show = label_symbols,
        non_label_count = adc_idx
    )
    l.xaxis[:rangeslider] = attr(;
        visible=slider, autorange=false, range=[0, dur(seq) * 1e3]
    )

    if !isempty(trigger_times)
        hover_positions = LinRange(0.0, 1.0, 51)
        trigger_x = repeat(trigger_times; inner=length(hover_positions))
        trigger_y = repeat(hover_positions, length(trigger_times))
        push!(p, scatter_fun(;
            x=trigger_x,
            y=trigger_y,
            name="Trigger",
            text=repeat(getproperty.(triggers, :label); inner=length(hover_positions)),
            hovertemplate="%{text}: %{x:.4f} ms<extra></extra>",
            xaxis=xaxis,
            yaxis="y3",
            showlegend=false,
            mode="markers",
            marker=attr(; color="rgba(255,0,0,0)", size=12),
        ))
        l.shapes = [
            attr(;
                type="line",
                x0=time,
                x1=time,
                y0=0,
                y1=1,
                xref=xaxis,
                yref="paper",
                layer="above",
                opacity=0.35,
                line=attr(; color="red", width=1),
            ) for time in trigger_times
        ]
        l.yaxis3 = attr(;
            overlaying="y",
            anchor=xaxis,
            range=[0.0, 1.0],
            fixedrange=true,
            visible=false,
        )
    end
    _add_physio!(p, l, scatter_fun, seq, physio, xaxis)
    return PlotlyBase.Plot(p, l; config)
end

function plot_M0_nonadaptive(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=false,
    show_seq_blocks=false,
    darkmode=false,
    range=[],
    title="",
)
    #Times
    seqd = KomaMRIBase.discretize(seq; sampling_rule=KomaMRIBase.MaxStepSizeRule(1, 5e-5))
    t, ts = seqd.t[1:(end - 1)], seqd.t[2:end]
    T0 = get_block_start_times(seq)
    #M0
    rf_idx, rf_types = KomaMRIBase.get_RF_types(seq, t)
    k, _ = KomaMRIBase.get_kspace(seqd; rf_idx, rf_types)
    #plots M0
    p = [scatter() for j in 1:4]
    p[1] = scatter(;
        x=ts * 1e3,
        y=k[:, 1],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)",
        name="M0x",
        legendgroup="Gx",
        marker=attr(; color="#636EFA"),
    )
    p[2] = scatter(;
        x=ts * 1e3,
        y=k[:, 2],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)",
        name="M0y",
        legendgroup="Gy",
        marker=attr(; color="#EF553B"),
    )
    p[3] = scatter(;
        x=ts * 1e3,
        y=k[:, 3],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)",
        name="M0z",
        legendgroup="Gz",
        marker=attr(; color="#00CC96"),
    )
    p[4] = scatter(;
        x=t[rf_idx] * 1e3,
        y=t[rf_idx] * 0,
        name="RF_center",
        marker=attr(; RF_CENTER_MARKER...),
        mode="markers",
        text=string.(rf_types)
    )
    #Layout and config
    l, config = generate_seq_time_layout_config(
        title, width, height, range, slider, show_seq_blocks, darkmode; T0
    )
    return PlotlyBase.Plot(p, l; config)
end

function plot_M1_nonadaptive(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=false,
    show_seq_blocks=false,
    darkmode=false,
    range=[],
    title="",
)
    #Times
    seqd = KomaMRIBase.discretize(seq; sampling_rule=KomaMRIBase.MaxStepSizeRule(1, 5e-5))
    t, ts = seqd.t[1:(end - 1)], seqd.t[2:end]
    T0 = get_block_start_times(seq)
    #M1
    rf_idx, rf_types = KomaMRIBase.get_RF_types(seq, t)
    k, _ = KomaMRIBase.get_M1(seqd; rf_idx, rf_types)
    #plots M1
    p = [scatter() for j in 1:4]
    p[1] = scatter(;
        x=ts * 1e3,
        y=k[:, 1],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms²)",
        name="M1x",
        legendgroup="Gx",
        marker=attr(; color="#636EFA"),
    )
    p[2] = scatter(;
        x=ts * 1e3,
        y=k[:, 2],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms²)",
        name="M1y",
        legendgroup="Gy",
        marker=attr(; color="#EF553B"),
    )
    p[3] = scatter(;
        x=ts * 1e3,
        y=k[:, 3],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms²)",
        name="M1z",
        legendgroup="Gz",
        marker=attr(; color="#00CC96"),
    )
    p[4] = scatter(;
        x=t[rf_idx] * 1e3,
        y=t[rf_idx] * 0,
        name="RF_center",
        marker=attr(; RF_CENTER_MARKER...),
        mode="markers",
        text=string.(rf_types)
    )
    #Layout and config
    l, config = generate_seq_time_layout_config(
        title, width, height, range, slider, show_seq_blocks, darkmode; T0
    )
    return PlotlyBase.Plot(p, l; config)
end

function plot_M2_nonadaptive(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=false,
    show_seq_blocks=false,
    darkmode=false,
    range=[],
    title="",
)
    #Times
    seqd = KomaMRIBase.discretize(seq; sampling_rule=KomaMRIBase.MaxStepSizeRule(1, 5e-5))
    t, ts = seqd.t[1:(end - 1)], seqd.t[2:end]
    T0 = get_block_start_times(seq)
    #M2
    rf_idx, rf_types = KomaMRIBase.get_RF_types(seq, t)
    k, _ = KomaMRIBase.get_M2(seqd; rf_idx, rf_types)
    #Plor M2
    p = [scatter() for j in 1:4]
    p[1] = scatter(;
        x=ts * 1e3,
        y=k[:, 1],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)",
        name="M2x",
        legendgroup="Gx",
        marker=attr(; color="#636EFA"),
    )
    p[2] = scatter(;
        x=ts * 1e3,
        y=k[:, 2],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)",
        name="M2y",
        legendgroup="Gy",
        marker=attr(; color="#EF553B"),
    )
    p[3] = scatter(;
        x=ts * 1e3,
        y=k[:, 3],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)",
        name="M2z",
        legendgroup="Gz",
        marker=attr(; color="#00CC96"),
    )
    p[4] = scatter(;
        x=t[rf_idx] * 1e3,
        y=t[rf_idx] * 0,
        name="RF_center",
        marker=attr(; RF_CENTER_MARKER...),
        mode="markers",
        text=string.(rf_types)
    )
    #Layout and config
    l, config = generate_seq_time_layout_config(
        title, width, height, range, slider, show_seq_blocks, darkmode; T0
    )
    return PlotlyBase.Plot(p, l; config)
end

function plot_eddy_currents_nonadaptive(
    seq::Sequence,
    λ;
    α=ones(size(λ)),
    width=nothing,
    height=nothing,
    slider=false,
    show_seq_blocks=false,
    darkmode=false,
    range=[],
    title="",
)
    #Times
    seqd = KomaMRIBase.discretize(seq + ADC(100, 100e-3); sampling_rule=KomaMRIBase.MaxStepSizeRule(1, 5e-5))
    t = seqd.t[2:end]
    T0 = get_block_start_times(seq)
    Gx, Gy, Gz = seqd.Gx[2:end], seqd.Gy[2:end], seqd.Gz[2:end]
    #Eddy currents per lambda
    Gec = zeros(length(t), 3)
    for (i, l) in enumerate(λ)
        aux, _ = KomaMRIBase.get_eddy_currents(seqd; λ=l)
        Gec .+= α[i] .* aux
    end
    #Plot eddy currents
    p = [scatter() for j in 1:4]
    p[1] = scatter(;
        x=t * 1e3,
        y=(Gx * 0 .+ Gec[:, 1]) * 1e3,
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
        name="ECx",
        legendgroup="Gx",
        marker=attr(; color="#636EFA"),
    )
    p[2] = scatter(;
        x=t * 1e3,
        y=(Gy * 0 .+ Gec[:, 2]) * 1e3,
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
        name="ECy",
        legendgroup="Gy",
        marker=attr(; color="#EF553B"),
    )
    p[3] = scatter(;
        x=t * 1e3,
        y=(Gz * 0 .+ Gec[:, 3]) * 1e3,
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
        name="ECz",
        legendgroup="Gz",
        marker=attr(; color="#00CC96"),
    )
    #Layout and config
    l, config = generate_seq_time_layout_config(
        title, width, height, range, slider, show_seq_blocks, darkmode; T0
    )
    return PlotlyBase.Plot(p, l; config)
end

function plot_slew_rate_nonadaptive(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=false,
    show_seq_blocks=false,
    darkmode=false,
    range=[],
    title="",
)
    #Times
    seqd = KomaMRIBase.discretize(seq; sampling_rule=KomaMRIBase.MaxStepSizeRule(1, 5e-5))
    ts = seqd.t
    T0 = get_block_start_times(seq)
    k, _ = KomaMRIBase.get_slew_rate(seqd)
    # Each slew value belongs to the interval ending at its timestamp.
    k = vcat(k[1:1, :], k)
    p = [scatter() for j in 1:4]
    p[1] = scatter(;
        x=ts * 1e3,
        y=k[:, 1],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m/ms)",
        name="SRx",
        legendgroup="Gx",
        mode="lines",
        line_shape="vh",
        marker=attr(; color="#636EFA"),
    )
    p[2] = scatter(;
        x=ts * 1e3,
        y=k[:, 2],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m/ms)",
        name="SRy",
        legendgroup="Gy",
        mode="lines",
        line_shape="vh",
        marker=attr(; color="#EF553B"),
    )
    p[3] = scatter(;
        x=ts * 1e3,
        y=k[:, 3],
        hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m/ms)",
        name="SRz",
        legendgroup="Gz",
        mode="lines",
        line_shape="vh",
        marker=attr(; color="#00CC96"),
    )
    #Layout and config
    l, config = generate_seq_time_layout_config(
        title, width, height, range, slider, show_seq_blocks, darkmode; T0
    )
    return PlotlyBase.Plot(p, l; config)
end


"""
    p = plot_image(image; height, width, zmin, zmax, darkmode, title)

Plot real image values directly. Volumes show linked x–y, z–y and x–z planes.
Complex images initially show magnitude and include
Magnitude/Phase buttons below the plot. Single-column images are shown as line profiles.

# Arguments
- `image`: a real or complex image matrix or volume

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `zmin`, `zmax`: color limits; for complex images these apply to magnitude only
- `colorscale`: (`="Greys"`) color scale for real values or magnitude;
  phase uses a cyclic scale over `[-π, π]` radians
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `title`: (`::String`, `=""`) plot title
- `spacing=(1,1,1)`: voxel spacing for volumes; axes display voxel indices

Middle-button drag adjusts window width horizontally and level vertically. Double-click resets window and zoom; phase keeps its fixed scale.

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the image matrix
"""
function plot_image(
    image;
    height=600,
    width=nothing,
    zmin=minimum(image[:]),
    zmax=maximum(image[:]),
    darkmode=false,
    title="",
    colorscale="Greys"
)
    #Layout
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
    l = Layout(;
        title=title,
        yaxis_title="y",
        xaxis_title="x",
        margin=attr(; t=50, l=0, r=0, b=0),
        yaxis=attr(; scaleanchor="x"),
        font_color=text_color,
        modebar=attr(;
            orientation="v", bgcolor=bgcolor, color=text_color, activecolor=plot_bgcolor
        ),
        xaxis=attr(; constrain="domain"),
        hovermode="closest",
        paper_bgcolor=bgcolor,
        plot_bgcolor=plot_bgcolor,
        xaxis_gridcolor=grid_color,
        yaxis_gridcolor=grid_color,
        xaxis_zerolinecolor=grid_color,
        yaxis_zerolinecolor=grid_color,
    )
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
    #Plot
    p = if size(image, 2) == 1
        l[:xaxis][:title] = "Index"
        l[:yaxis] = attr(; title="Magnitude", range=[zmin, zmax], scaleanchor=false,
            gridcolor=grid_color, zerolinecolor=grid_color)
        l[:margin] = attr(; t=50, l=60, r=20, b=50)
        scatter(; x=0:(length(image) - 1), y=vec(image), mode="lines", line_color="#2a7fb8")
    else
        heatmap(; z=collect(eachcol(image)), transpose=false, zmin=zmin, zmax=zmax,
            meta=Dict(:image_window=>[zmin,zmax]), colorscale=colorscale, hoverinfo="x+y+z")
    end
    config = PlotConfig(;
        displaylogo=false,
        toImageButtonOptions=attr(;
            format="svg", # one of png, svg, jpeg, webp
        ).fields,
        modeBarButtonsToRemove=[
            "zoom",
            "autoScale",
            "resetScale2d",
            "pan",
            "tableRotation",
            "resetCameraLastSave",
            "zoomIn",
            "zoomOut",
        ],
    )
    return pluto_plot(PlotlyBase.Plot(p, l; config))
end

function plot_image(image::AbstractArray{<:Complex}; kwargs...)
    magnitude = plot_image(abs.(image); kwargs...)
    phase_limit = Float64(π)
    phase = plot_image(angle.(image); kwargs..., zmin=-phase_limit, zmax=phase_limit, colorscale=PHASE_COLORSCALE)
    profile = ndims(image) < 3 && size(image, 2) == 1
    if profile
        lower, upper = magnitude.layout[:yaxis][:range]
        magnitude.layout[:yaxis][:range] = [lower, 1.1 * upper]
        phase.layout[:yaxis][:title] = "Phase (rad)"
    else
        for trace in phase.data
            delete!(trace.fields, :meta)
            colorbar = get!(trace.fields,:colorbar,attr())
            colorbar[:title] = attr(; text="rad")
            colorbar[:tickvals] = [-phase_limit,-phase_limit/2,0,phase_limit/2,phase_limit]
            colorbar[:ticktext] = ["−π","−π/2","0","π/2","π"]
        end
    end
    components = (magnitude, phase)
    traces = reduce(vcat, [component.data for component in components])
    count = length(magnitude.data)
    for (index, trace) in enumerate(traces)
        trace[:visible] = index <= count
    end
    layout = magnitude.layout
    if haskey(get(layout,:meta,Dict()),"image_slices")
        model = layout[:meta]["image_slices"]
        layout[:meta]["image_slices"] = (; model...,values=[model.values;phase.layout[:meta]["image_slices"].values])
        for (index,trace) in enumerate(traces)
            trace[:uid] = "image-component-$index"
        end
    end
    layout[:showlegend] = false
    layout[:margin][:b] = 100
    layout[:updatemenus] = [attr(;
        name="image-component", type="buttons", direction="right",
        x=0.5, xanchor="center", y=0, yanchor="top", pad=attr(; t=50),
        bgcolor="#2a7fb8", bordercolor="#2a7fb8", font=attr(; color="white"),
        showactive=true, active=0,
        buttons=[attr(; label, method="update", args=[
            attr(; visible=[fill(index == 1,count); fill(index == 2,count)]),
            profile ? Dict("yaxis.range" => components[index].layout[:yaxis][:range],
                "yaxis.title" => components[index].layout[:yaxis][:title]) : attr(),
        ]) for (index, label) in enumerate(("Magnitude", "Phase"))],
    )]
    return pluto_plot(PlotlyBase.Plot(traces, layout; config=magnitude.config))
end

"""
    p = plot_kspace(seq::Sequence; width=nothing, height=nothing, darkmode=false, view_2d=false, adaptive=false, adaptive_navigator=false, physio=NoPhysioSignal())

Plots the k-space of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `view_2d`: (`::Bool`, `=false`) boolean to indicate whether to use a 2D kx/ky plot
- `physio=NoPhysioSignal()`: physiological signal used to resolve trigger waits and ADC times
- `adaptive=false`: use an experimental live `SpatialPlot` when true. Julia retains
  the full trajectory; zoom/rotation fetch bounded sets of original segments and ADCs.
- `adaptive_navigator=false`: show a compact adaptive sequence navigator when true, without
  a legend or secondary slider. Zoom its time axis to filter k-space after dragging
  stops; colors retain their global ADC sample count. The sequence navigator requires a running Julia
  session in both k-space modes.

# Returns
- `p`: `PlotlyBase.Plot` with `adaptive=false`, or live `SpatialPlot` otherwise.
  Adaptive updates require the Julia session to remain running.

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_kspace(seq)
```
"""
function plot_kspace(seq::Sequence; width=nothing, height=nothing, darkmode=false, view_2d=false,
    adaptive=false, adaptive_navigator=false, physio=NoPhysioSignal())
    seq = resolve_triggers(seq, physio)
    (adaptive || adaptive_navigator) && return kspace_spatial_plot(seq; width, height, darkmode, view_2d, adaptive, adaptive_navigator)
    #Calculations of theoretical k-space
    kspace, kspace_adc = get_kspace(seq; sampling_rule=KomaMRIBase.MaxStepSizeRule(1, 5e-5))
    t_adc = KomaMRIBase.get_adc_sampling_times(seq)
    return pluto_plot(kspace_plot(kspace, kspace_adc, t_adc, length(t_adc); width, height, darkmode, view_2d))
end

function kspace_plot(kspace, kspace_adc, t_adc, num_adc; width, height, darkmode, view_2d)
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
    c = "gray"
    #Layout
    mink = minimum(kspace_adc; dims=1)
    maxk = maximum(kspace_adc; dims=1)
    dW = maximum(maxk .- mink; dims=2) * 0.3
    mink .-= dW
    maxk .+= dW
    modebar = attr(;
        orientation="h",
        yanchor="bottom",
        xanchor="right",
        y=1,
        x=0,
        bgcolor=bgcolor,
        color=text_color,
        activecolor=plot_bgcolor,
    )
    legend = attr(; orientation="h", yref="container", yanchor="top", xanchor="left", y=1, x=0)
    if view_2d
        l = Layout(;
            paper_bgcolor=bgcolor,
            plot_bgcolor=plot_bgcolor,
            xaxis=attr(;
                title="kx [m⁻¹]",
                range=[mink[1], maxk[1]],
                gridcolor=grid_color,
                zerolinecolor=grid_color,
                scaleanchor="y",
            ),
            yaxis=attr(;
                title="ky [m⁻¹]",
                range=[mink[2], maxk[2]],
                gridcolor=grid_color,
                zerolinecolor=grid_color,
                scaleratio=1,
            ),
            modebar,
            legend,
            font_color=text_color,
            margin=attr(; t=0, l=0, r=0),
        )
        p = [
            scattergl(;
                x=kspace[:, 1],
                y=kspace[:, 2],
                mode="lines",
                line=attr(; color=c),
                name="Trajectory",
                hoverinfo="skip",
            ),
            scattergl(;
                x=kspace_adc[:, 1],
                y=kspace_adc[:, 2],
                text=round.(t_adc * 1e3, digits=3),
                mode="markers",
                marker=attr(; size=4),
                name="ADC",
                hovertemplate="ADC count: %{marker.color}<br>kx: %{x:.1f} m⁻¹<br>ky: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{text} ms<extra></extra>",
            ),
            scattergl(;
                x=[0], y=[0], mode="markers", name="k=0",
                marker=attr(; symbol="cross", size=10, color="red"),
            ),
        ]
        modebar_buttons = ["zoom", "pan", "resetScale2d", "autoScale2d"]
    else
        l = Layout(;
            paper_bgcolor=bgcolor,
            scene=attr(;
                xaxis=attr(;
                    title="kx [m⁻¹]",
                    range=[mink[1], maxk[1]],
                    backgroundcolor=plot_bgcolor,
                    gridcolor=grid_color,
                    zerolinecolor=grid_color,
                ),
                yaxis=attr(;
                    title="ky [m⁻¹]",
                    range=[mink[2], maxk[2]],
                    backgroundcolor=plot_bgcolor,
                    gridcolor=grid_color,
                    zerolinecolor=grid_color,
                ),
                zaxis=attr(;
                    title="kz [m⁻¹]",
                    range=[mink[3], maxk[3]],
                    backgroundcolor=plot_bgcolor,
                    gridcolor=grid_color,
                    zerolinecolor=grid_color,
                ),
            ),
            modebar,
            legend,
            font_color=text_color,
            scene_camera_eye=attr(; x=0, y=0, z=1.7),
            scene_camera_up=attr(; x=0, y=1.0, z=0),
            scene_aspectmode="cube",
            margin=attr(; t=0, l=0, r=0),
        )
        p = [
            scatter3d(;
                x=kspace[:, 1],
                y=kspace[:, 2],
                z=kspace[:, 3],
                mode="lines",
                line=attr(; color=c),
                name="Trajectory",
                hoverinfo="skip",
            ),
            scatter3d(;
                x=kspace_adc[:, 1],
                y=kspace_adc[:, 2],
                z=kspace_adc[:, 3],
                text=round.(t_adc * 1e3, digits=3),
                mode="markers",
                marker=attr(; size=2),
                name="ADC",
                hovertemplate="ADC count: %{marker.color}<br>kx: %{x:.1f} m⁻¹<br>ky: %{y:.1f} m⁻¹<br>kz: %{z:.1f} m⁻¹<br><b>t_acq</b>: %{text} ms<extra></extra>",
            ),
            scatter3d(;
                x=[0], y=[0], z=[0], mode="markers", name="k=0",
                marker=attr(; symbol="cross", size=10, color="red"),
            ),
        ]
        modebar_buttons = [
            "zoom",
            "pan",
            "tableRotation",
            "resetCameraLastSave3d",
            "orbitRotation",
            "resetCameraDefault3d",
        ]
    end
    p[2][:marker][:color] = collect(0:length(t_adc)-1)
    p[2][:marker][:colorscale] = KSPACE_COLORSCALE
    p[2][:marker][:cmin] = 0
    p[2][:marker][:cmax] = max(num_adc - 1, 1)
    p[2][:marker][:showscale] = true
    p[2][:marker][:colorbar] = attr(; title=attr(; text="ADC count", side="right"), thickness=12, len=0.6,
        tickformat="~s", outlinewidth=0)
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
    config = PlotConfig(;
        displaylogo=false,
        toImageButtonOptions=attr(;
            format="svg", # one of png, svg, jpeg, webp
        ).fields,
        modeBarButtonsToRemove=modebar_buttons,
    )
    return PlotlyBase.Plot(p, l; config)
end

"""
    p = plot_phantom_map(obj::Phantom, key::Symbol; kwargs...)

Plots a phantom map for a specific spin parameter given by `key`.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `key`: (`::Symbol`, opts: [`:ρ`, `:T1`, `:T2`, `:T2s`, `:x`, `:y`, `:z`, `:speed`, `:vx`, `:vy`, `:vz`]) symbol for
    displaying different parameters of the phantom spins
    (`:speed` and signed velocity components use frame-to-frame displacement in cm/s;
    components use blue/white/red for negative/zero/positive with fixed symmetric limits.
    Exiting spins flash magenta before a reset and cyan upon re-entry;
    reset jumps are excluded from the velocity scale).

# Keywords — both modes
- `adaptive=false`: return a self-contained `PlotlyBase.Plot`; `true` returns a live
  [`SpatialPlot`](@ref), requiring Julia for zoom, camera and time updates.
- `height`: (`::Integer`, `=700`) plot height
- `width`: (`::Integer`, `=nothing`) plot width
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `view_2d`: use a 2D scatter plot; defaults to true for fewer than three spatial dimensions
- `colorbar`: (`::Bool`, `=true`) boolean to indicate whether to display a colorbar
- `time_samples=0`: add an evenly spaced time grid to the native motion timepoints;
  duplicates are removed, and native timepoints are never discarded. Static phantoms
  have one frame. Nonadaptive plots embed every resulting frame.

# Keywords — only `adaptive=false`
- `max_spins`: cap the embedded spin subset; defaults to 20,000 when omitted.
  Passing it with `adaptive=true` throws `ArgumentError`.

# Returns
- `p`: `PlotlyBase.Plot` when `adaptive=false`, otherwise [`SpatialPlot`](@ref)

# References
Colormaps from https://github.com/markgriswold/MRFColormaps
Towards Unified Colormaps for Quantitative MRF Data, Mark Griswold, et al. (2018).

# Examples
```julia-repl
julia> obj2D, obj3D = brain_phantom2D(), brain_phantom3D();

julia> plot_phantom_map(obj2D, :ρ)

julia> plot_phantom_map(obj3D, :ρ)
```
"""
function plot_phantom_map(obj::Phantom, key::Symbol; height=700, width=nothing,
    darkmode=false, view_2d=sum(KomaMRIBase.get_dims(obj)) < 3, colorbar=true,
    max_spins=nothing, adaptive=false,
    time_samples=0, kwargs...)
    return phantom_plot(obj, key; adaptive, max_spins, time_samples,
        height, width, darkmode, view_2d, colorbar, kwargs...)
end



# Nonadaptive signal figure, shared with the adaptive signal source.
function plot_signal_nonadaptive(
    raw::RawAcquisitionData;
    width=nothing,
    height=nothing,
    slider=false,
    show_sim_blocks=false,
    darkmode=false,
    range=[],
    gl=false,
)
    not_Koma = raw.params["systemVendor"] != "KomaMRI.jl"
    t = []
    ncoils = size(first(raw.profiles).data, 2)
    scatter_fun = gl ? scattergl : scatter
    components = ((abs, "|S(t)|", scatter), (real, "Re{S(t)}", scatter_fun),
        (imag, "Im{S(t)}", scatter_fun))
    trace_coils = repeat(1:ncoils; inner=length(components))
    signal_type = Union{Missing,eltype(first(raw.profiles).data)}
    signals = [signal_type[] for _ in 1:ncoils]
    current_t0 = 0
    for p in raw.profiles
        dt = p.head.sample_time_us != 0 ? p.head.sample_time_us * 1e-3 : 1
        t0 = p.head.acquisition_time_stamp * 1e-3 #This parameter is used in Koma to store the time offset
        N = p.head.number_of_samples != 0 ? p.head.number_of_samples : 1
        if not_Koma
            t0 = current_t0 * dt
            current_t0 += N
        end
        if N != 1
            append!(t, t0 .+ (0:dt:(dt * (N - 1))))
        else
            append!(t, t0)
        end
        for coil in 1:ncoils
            append!(signals[coil], @view p.data[:, coil])
            push!(signals[coil], missing)
        end
        #To generate gap
        append!(t, t[end])
    end
    ymin, ymax = extrema(
        value for signal in signals for sample in skipmissing(signal) for
        value in (real(sample), imag(sample), abs(sample))
    )
    padding = iszero(ymax - ymin) ? max(abs(ymax), one(ymax)) : (ymax - ymin) / 20
    signal_range = [ymin - padding, ymax + padding]
    #Show simulation blocks
    shapes = []
    annotations = []
    type_names = ["precession", "excitation"]
    if !not_Koma && show_sim_blocks
        t_sim_parts = raw.params["userParameters"]["t_sim_parts"]
        type_sim_parts = raw.params["userParameters"]["type_sim_parts"]

        current_type = -1
        for i in eachindex(t_sim_parts[1:(end - 1)])
            aux = rect(;
                xref="x",
                yref="paper",
                x0=t_sim_parts[i] * 1e3,
                y0=0,
                x1=t_sim_parts[i + 1] * 1e3,
                y1=1,
                fillcolor=type_sim_parts[i] ? "Purple" : "Blue",
                opacity=0.1,
                layer="below",
                line_width=2,
            )
            push!(shapes, aux)

            if type_sim_parts[i] != current_type
                aux = attr(;
                    xref="x",
                    yref="paper",
                    x=t_sim_parts[i] * 1e3,
                    y=1,
                    showarrow=false,
                    text=type_names[type_sim_parts[i] + 1],
                )
                push!(annotations, aux)
                current_type = type_sim_parts[i]
            end
        end
    end
    #PLOT
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
    l = Layout(;
        hovermode="closest",
        xaxis_title="",
        modebar=attr(;
            orientation="h",
            yanchor="bottom",
            xanchor="right",
            y=1,
            x=0,
            bgcolor=bgcolor,
            color=text_color,
            activecolor=plot_bgcolor,
        ),
        legend=attr(; orientation="h", yref="container", yanchor="top", xanchor="left", y=1, x=0),
        plot_bgcolor=plot_bgcolor,
        paper_bgcolor=bgcolor,
        xaxis_gridcolor=grid_color,
        yaxis_gridcolor=grid_color,
        xaxis_zerolinecolor=grid_color,
        yaxis_zerolinecolor=grid_color,
        yaxis_range=signal_range,
        font_color=text_color,
        yaxis_fixedrange=false,
        yaxis_automargin=false,
        xaxis=attr(;
            automargin=false,
            ticksuffix=" ms",
            range=range[:],
            rangeslider=attr(; visible=slider),
            rangeselector=attr(;
                buttons=[
                    attr(; count=1, label="1m", step=10, stepmode="backward"),
                    attr(; step="all"),
                ],
            ),
        ),
        shapes=shapes,
        annotations=annotations,
        sliders=ncoils == 1 ? [] : [
            attr(;
                active=0,
                currentvalue=attr(; prefix="Coil: "),
                steps=[
                    attr(;
                        label=string(coil),
                        method="restyle",
                        args=[attr(; visible=trace_coils .== coil)],
                    ) for coil in 1:ncoils
                ],
                x=0.1,
                len=0.8,
                y=1.12,
                yanchor="bottom",
            ),
        ],
        margin=attr(; t=ncoils == 1 ? 6 : 125, l=30, r=6, b=24),
    )
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
    p = [
        plot_component(;
            x=t,
            y=component.(signal),
            name,
            hovertemplate="(%{x:.4f} ms, %{y:.3f} a.u.)",
            visible=coil == 1,
        ) for (coil, signal) in enumerate(signals)
        for (component, name, plot_component) in components
    ]
    config = PlotConfig(;
        displaylogo=false,
        toImageButtonOptions=attr(;
            format="svg", # one of png, svg, jpeg, webp
        ).fields,
        modeBarButtonsToRemove=[
            "zoom",
            "autoScale",
            "resetScale2d",
            "pan",
            "tableRotation",
            "resetCameraLastSave",
            "zoomIn",
            "zoomOut",
        ],
        # modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
    )
    return PlotlyBase.Plot(p, l; config)
end


"""
    str = plot_dict(dict::Dict)

Generates an HTML table based on the dictionary `dict`.

# Arguments
- `dict`: (`::Dict`) dictionary

# Returns
- `str`: (`::String`) dictionary as an HTML table
"""
function plot_dict(dict::Dict)
    html = """
    <table class="table table-dark table-striped">
        <thead>
            <tr>
            <th scope="col">#</th>
            <th scope="col">Name</th>
            <th scope="col">Value</th>
            </tr>
        </thead>
        <tbody>
    """
    i = 1
    for (key, val) in dict
        html *= """
            <tr>
                <th scope="row">$i</th>
                <td>$(string(key))</td>
                <td>$(string(val))</td>
            </tr>
        """
        i += 1
    end
    return html *= "</tbody></table>"
end

function plot_seqd_marker_symbols(seq, seqd, sampling_rule; freq_in_phase=false)
    hasproperty(sampling_rule, :rule) || return fill(:circle, length(seqd.t))
    boundary_seqd = KomaMRIBase.discretize(seq; sampling_rule=getproperty(sampling_rule, :rule), freq_in_phase)
    boundary_t = Set(boundary_seqd.t)
    return [t in boundary_t ? :circle : Symbol("line-ns") for t in seqd.t]
end

# Nonadaptive sampled-waveform figure, shared with the adaptive source.
function plot_seqd_nonadaptive(seq::Sequence; sampling_rule=KomaMRIBase.MaxStepSizeRule(1e-3, 5e-5), show_rf_frame=true, freq_in_phase=false)
    seqd = KomaMRIBase.discretize(seq; sampling_rule, freq_in_phase)
    marker_symbol = plot_seqd_marker_symbols(seq, seqd, sampling_rule; freq_in_phase)
    marker_line_width = [s == :circle ? 0 : 2 for s in marker_symbol]
    is_real_rf = all(x -> isapprox(imag(x), 0; atol=eps()), seqd.B1)
    B1 = is_real_rf ? real.(seqd.B1) : abs.(seqd.B1)
    B1_phase = is_real_rf ? zero.(real.(seqd.B1)) : angle.(seqd.B1)
    Gx = scattergl(;
        x=seqd.t * 1e3,
        y=seqd.Gx * 1e3,
        name="Gx",
        mode="markers+lines",
        marker_symbol,
        legendgroup="Gx",
        marker=attr(; color="#636EFA", size=8, line=attr(; color="#636EFA", width=marker_line_width)),
        line=attr(; color="#636EFA"),
    )
    Gy = scattergl(;
        x=seqd.t * 1e3,
        y=seqd.Gy * 1e3,
        name="Gy",
        mode="markers+lines",
        marker_symbol,
        legendgroup="Gy",
        marker=attr(; color="#EF553B", size=8, line=attr(; color="#EF553B", width=marker_line_width)),
        line=attr(; color="#EF553B"),
    )
    Gz = scattergl(;
        x=seqd.t * 1e3,
        y=seqd.Gz * 1e3,
        name="Gz",
        mode="markers+lines",
        marker_symbol,
        legendgroup="Gz",
        marker=attr(; color="#00CC96", size=8, line=attr(; color="#00CC96", width=marker_line_width)),
        line=attr(; color="#00CC96"),
    )
    B1_abs = scattergl(;
        x=seqd.t * 1e3,
        y=B1 * 1e6,
        name="|B1|_AM",
        mode="markers+lines",
        marker_symbol,
        legendgroup="|B1|_AM",
        marker=attr(; color="#AB63FA", size=8, line=attr(; color="#AB63FA", width=marker_line_width)),
        line=attr(; color="#AB63FA"),
    )
    B1_angle = scattergl(;
        x=seqd.t * 1e3,
        y=B1_phase,
        name="∠B1_AM",
        mode="markers+lines",
        marker_symbol,
        legendgroup="∠B1_AM",
        marker=attr(; color="#FFA15A", size=8, line=attr(; color="#FFA15A", width=marker_line_width)),
        line=attr(; color="#FFA15A"),
    )
    ADC = scattergl(;
        x=seqd.t[seqd.ADC] * 1e3,
        y=zeros(sum(seqd.ADC)),
        name="ADC",
        mode="markers",
        marker_symbol=:x,
        legendgroup="ADC",
        marker=attr(; color="#19D3F3"),
    )
    B1_Δf = scattergl(;
        x=seqd.t * 1e3,
        y=seqd.Δf * 1e-3,
        name="Δf_FM",
        mode="markers+lines",
        marker_symbol,
        visible="legendonly",
        legendgroup="Δf_FM",
        marker=attr(; color="#AB63FA", size=8, line=attr(; color="#AB63FA", width=marker_line_width)),
        line=attr(; color="#AB63FA"),
    )
    excitation_bool = scattergl(;
        x=seqd.t * 1e3,
        y=Float64.([seqd.excitation_bool; false] .| [false; seqd.excitation_bool]),
        name="excitation_bool",
        mode="lines",
        visible="legendonly",
        legendgroup="excitation_bool",
        line=attr(; color="#AB63FA", dash="dot"),
    )
    p = [Gx, Gy, Gz, B1_abs, B1_angle, ADC, B1_Δf, excitation_bool]
    if show_rf_frame
        push!(p, scattergl(;
            x=seqd.t * 1e3,
            y=seqd.ψ,
            name="ψ_FM",
            mode="markers+lines",
            marker_symbol,
            legendgroup="ψ_FM",
            marker=attr(; color="#FF6692", size=8, line=attr(; color="#FF6692", width=marker_line_width)),
            line=attr(; color="#FF6692"),
            visible="legendonly",
        ))
    end
    return PlotlyBase.Plot(p)
end

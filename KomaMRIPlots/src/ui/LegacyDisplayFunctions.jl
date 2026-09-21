const RF_CENTER_MARKER = (color="#FF0000", symbol="x")

"""
    p = plot_seq_legacy(seq::Sequence; kwargs...)

Plots a sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title
- `freq_in_phase`: (`::Bool`, `=true`) Include FM modulation in RF phase
- `show_rf_frame`: (`::Bool`, `=false`) plot RF rotating-frame phase
- `gl`: (`::Bool`, `=false`) use the Plotly `scattergl` trace (faster)
- `max_rf_samples`: (`::Integer`, `=100`) maximum number of RF samples
- `show_adc`: (`::Bool`, `=false`) plot ADC samples with markers
- `physio`: (`=NoPhysioSignal()`) physiological signal used to resolve triggers

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_seq_legacy(seq)
```
"""
function plot_seq_legacy(
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
    trigger_times = [
        (T0[i] + delay) * 1e3
        for (i, extensions) in enumerate(seq.EXT)
        for extension in extensions
        for delay in _trigger_delays(extension)
    ]
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
            hovertemplate="Trigger: %{x:.4f} ms<extra></extra>",
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

"""
    p = plot_M0_legacy(seq::Sequence; kwargs...)

Plots the zero order moment (M0) of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the moment M0 of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_M0_legacy(seq)
```
"""
function plot_M0_legacy(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=true,
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

"""
    p = plot_M1_legacy(seq::Sequence; kwargs...)

Plots the first order moment (M1) of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the moment M1 of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_M1_legacy(seq)
```
"""
function plot_M1_legacy(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=true,
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

"""
    p = plot_M2_legacy(seq::Sequence; kwargs...)

Plots the second order moment (M2) of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the moment M2 of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_M2_legacy(seq)
```
"""
function plot_M2_legacy(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=true,
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

"""
    p = plot_eddy_currents_legacy(seq::Sequence, λ; kwargs...)

Plots the eddy currents of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `λ`: (`::Real`, `[s]`) time constant for the decay of Eddy currents

# Keywords
- `α`: (`::Vector{Real}`, `=ones(size(λ))`) eddy currents factors
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the Eddy currents of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_eddy_currents_legacy(seq, 80e-3)
```
"""
function plot_eddy_currents_legacy(
    seq::Sequence,
    λ;
    α=ones(size(λ)),
    width=nothing,
    height=nothing,
    slider=true,
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

"""
    p = plot_slew_rate_legacy(seq::Sequence; kwargs...)

Plots the slew rate currents of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the slew rate currents of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_slew_rate_legacy(seq)
```
"""
function plot_slew_rate_legacy(
    seq::Sequence;
    width=nothing,
    height=nothing,
    slider=true,
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
    p = plot_signal_legacy(raw::RawAcquisitionData; kwargs...)

Plots a raw signal in ISMRMRD format.

For multi-coil data, a slider selects the displayed receive channel.

# Arguments
- `raw`: (`::RawAcquisitionData`) RawAcquisitionData struct (raw signal in ISMRMRD format)

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_sim_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the raw signal

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq");

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> raw = simulate(obj, seq, sys)

julia> plot_signal_legacy(raw)
```
"""
function plot_signal_legacy(
    raw::RawAcquisitionData;
    width=nothing,
    height=nothing,
    slider=true,
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
        legend=attr(; orientation="h", yanchor="bottom", xanchor="left", y=1, x=0),
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
    p = plot_seqd_legacy(seq::Sequence; sampling_rule=KomaMRIBase.MaxStepSizeRule(1e-3, 5e-5))

Plots a sampled sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `sampling_rule`: controls how the sequence sampling grid is refined
- `show_rf_frame`: (`::Bool`, `=true`) plot RF rotating-frame phase
- `freq_in_phase`: (`::Bool`, `=false`) fold RF frequency modulation into the complex RF waveform

# Returns
- `p`: (`::PlotlyBase.Plot`) plot of the sampled Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_seqd_legacy(seq)
```
"""
function plot_seqd_legacy(seq::Sequence; sampling_rule=KomaMRIBase.MaxStepSizeRule(1e-3, 5e-5), show_rf_frame=true, freq_in_phase=false)
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

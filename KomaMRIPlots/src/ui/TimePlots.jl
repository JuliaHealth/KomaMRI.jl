struct TracePlotSource{D,L,C,I,R}
    data::D
    layout::L
    config::C
    samples::I
    full_range::R
end

function TracePlotSource(plot)
    samples = [index_samples(get(t, :x, Float64[]), get(t, :y, Float64[])) for t in plot.data]
    times = [t[:x][i] for (t, s) in zip(plot.data, samples) for i in s.valid]
    limits = isempty(times) ? [0.0, 1.0] : collect(extrema(times))
    limits[1] == limits[2] && (limits .+= [-0.5, 0.5])
    return TracePlotSource(plot.data, plot.layout, plot.config, samples, limits)
end

plot_template(source::TracePlotSource) = source

sample_field(value, indices, n) = value
sample_field(value::AbstractVector, indices, n) = length(value) == n ?
    [i == 0 ? nothing : value[i] for i in indices] : value
sample_field(value::AbstractDict, indices, n) =
    Dict(key => sample_field(item, indices, n) for (key, item) in pairs(value))
sample_field(value::PlotlyBase.AbstractPlotlyAttribute, indices, n) = sample_field(value.fields, indices, n)

function window_data(source::TracePlotSource, interval, width, visibility)
    bins = clamp(floor(Int, width - 70), 24, 1000)
    shown = available = 0
    data = map(eachindex(source.data)) do j
        trace = source.data[j]
        x, y = get(trace, :x, Float64[]), get(trace, :y, Float64[])
        indices, count = visibility[j] === true ?
            selected_samples(x, y, source.samples[j], interval, bins) : (Int[], 0)
        available += count
        shown += Base.count(i -> i != 0 && interval[1] <= x[i] <= interval[2], indices)
        visibility[j] == "legendonly" && (indices = [0])
        fields = copy(trace.fields)
        for name in (:x, :y, :text, :hovertext, :customdata, :ids, :marker)
            haskey(fields, name) && (fields[name] = sample_field(fields[name], indices, length(x)))
        end
        fields[:mode] = get(fields, :mode, "lines")
        fields[:visible] = visibility[j]
        fields[:uid] = "trace-$j"
        fields
    end
    layout = deepcopy(source.layout)
    if haskey(layout, :shapes)
        layout[:shapes] = filter(get(layout, :shapes, [])) do shape
            get(shape, :xref, "x") != "x" ||
                (get(shape, :x1, interval[2]) >= interval[1] && get(shape, :x0, interval[1]) <= interval[2])
        end
    end
    return (; data, layout, source.config, source.full_range, shown, available)
end

"""
    TimePlot

A Plotly time-series viewer that obtains visible samples from Julia on zoom or pan.
Display it directly in Pluto, IJulia/Jupyter, or a browser, or embed it in a Bonito app.
The Julia session must remain running for new detail. Created with `adaptive=true`;
the default `adaptive=false` returns a regular `PlotlyBase.Plot`.
Keyword availability is documented on each `plot_*` function; nonadaptive-only
keywords cannot be passed with `adaptive=true`.
"""
struct TimePlot{S}
    source::S
end

"""
    plot_seq(seq; adaptive=false, physio=NoPhysioSignal(), slider=false, kwargs...)

Plot sequence waveforms as a regular `PlotlyBase.Plot`. Set `adaptive=true` for a
live [`TimePlot`](@ref) that fetches detail from Julia on zoom or pan. The Julia
session must remain running; an indicator marks reduced samples.

# Keywords — both modes
- `width=nothing`, `height=nothing`: plot dimensions in pixels.
- `slider=false`, `darkmode=false`, `title=""`: display controls.
- `range=[]`: initial time range in milliseconds.
- `physio=NoPhysioSignal()`: resolve triggers and display physiology.

# Keywords — only `adaptive=false`
- `gl=false`: use WebGL traces.
- `show_seq_blocks=false`: label sequence block boundaries.
- `show_adc=false`: show individual ADC sample markers.
- `max_rf_samples=100`: cap embedded RF samples per event.
- `freq_in_phase=false`, `show_rf_frame=false`: RF phase display.
- `xaxis="x"`, `yaxis="y"`, `showlegend=true`: trace layout controls.

Passing nonadaptive-only keywords with `adaptive=true` throws `ArgumentError`.
Adaptive RF samples are uncapped; vertical ADC ticks appear automatically when
sufficiently separated. Trigger timing and labels use the full resolved sequence.
Multiple RF channels currently use full trace construction before adaptive sampling.
"""
function plot_seq(seq::Sequence; physio=NoPhysioSignal(), width=nothing, height=nothing,
    slider=false, darkmode=false, range=[], title="", adaptive=false, kwargs...)
    adaptive || return plot_seq_nonadaptive(seq; physio, width, height, slider, darkmode,
        range, title, kwargs...)
    isempty(kwargs) || throw(ArgumentError("These keywords require adaptive=false: $(join(keys(kwargs), ", "))."))
    source = if size(seq.RF, 1) == 1
        sequence_source(seq; physio, width, height, slider, darkmode, range, title)
    else
        TracePlotSource(plot_seq_nonadaptive(seq; physio, width, height, slider, darkmode,
            range, title, max_rf_samples=typemax(Int)))
    end
    return TimePlot(source)
end

time_plot(plot, adaptive) = adaptive ? TimePlot(TracePlotSource(plot)) : plot

"""
    plot_M0(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...)

Plot the zero-order gradient moment. Return a self-contained `PlotlyBase.Plot`,
or a live [`TimePlot`](@ref) with `adaptive=true`.

# Keywords — both modes
`width=nothing`, `height=nothing`, `slider=true`, `show_seq_blocks=false`,
`darkmode=false`, `range=[]` (milliseconds), `title=""`, and `physio=NoPhysioSignal()`.
Both modes construct the full moment arrays; adaptive mode limits displayed samples.
There are no nonadaptive-only keywords.
"""
plot_M0(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...) =
    time_plot(plot_M0_nonadaptive(resolve_triggers(seq, physio); kwargs...), adaptive)
"""
    plot_M1(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...)

Plot the first-order gradient moment. Keywords and mode behavior match [`plot_M0`](@ref).
All plotting keywords apply in both modes.
"""
plot_M1(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...) =
    time_plot(plot_M1_nonadaptive(resolve_triggers(seq, physio); kwargs...), adaptive)
"""
    plot_M2(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...)

Plot the second-order gradient moment. Keywords and mode behavior match [`plot_M0`](@ref).
All plotting keywords apply in both modes.
"""
plot_M2(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...) =
    time_plot(plot_M2_nonadaptive(resolve_triggers(seq, physio); kwargs...), adaptive)
"""
    plot_slew_rate(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...)

Plot gradient slew rate. Keywords and mode behavior match [`plot_M0`](@ref).
All plotting keywords apply in both modes.
"""
plot_slew_rate(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...) =
    time_plot(plot_slew_rate_nonadaptive(resolve_triggers(seq, physio); kwargs...), adaptive)
"""
    plot_eddy_currents(seq, λ; adaptive=false, physio=NoPhysioSignal(), kwargs...)

Plot eddy currents with time constants `λ` and weights `α=ones(size(λ))`.
Other keywords and mode behavior match [`plot_M0`](@ref). All plotting keywords,
including `α`, apply in both modes.
"""
plot_eddy_currents(seq, λ; adaptive=false, physio=NoPhysioSignal(), kwargs...) =
    time_plot(plot_eddy_currents_nonadaptive(resolve_triggers(seq, physio), λ; kwargs...), adaptive)
"""
    plot_seqd(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...)

Plot discretized waveforms as a self-contained `PlotlyBase.Plot`, or a live
[`TimePlot`](@ref) with `adaptive=true`.

# Keywords — both modes
`sampling_rule=MaxStepSizeRule(1e-3, 5e-5)`, `show_rf_frame=true`,
`freq_in_phase=false`, and `physio=NoPhysioSignal()`.
`sampling_rule` controls the underlying discretization, not viewport downsampling.
Both modes construct the full sampled arrays. There are no nonadaptive-only keywords.
"""
plot_seqd(seq; adaptive=false, physio=NoPhysioSignal(), kwargs...) =
    time_plot(plot_seqd_nonadaptive(resolve_triggers(seq, physio); kwargs...), adaptive)
"""
    plot_signal(raw; adaptive=false, kwargs...)

Plot raw signals with coil selection as a self-contained `PlotlyBase.Plot`, or a live
[`TimePlot`](@ref) with `adaptive=true`.

# Keywords — both modes
`width=nothing`, `height=nothing`, `slider=true`, `show_sim_blocks=false`,
`darkmode=false`, `range=[]` (milliseconds), and `gl=false`.
Both modes construct the full signal arrays. There are no nonadaptive-only keywords.
"""
plot_signal(raw; adaptive=false, kwargs...) = time_plot(plot_signal_nonadaptive(raw; kwargs...), adaptive)

const TIME_PLOTLY = Bonito.Asset(joinpath(artifact"plotly-artifacts", "plotly.min.js"); name="Plotly")
const TIME_PLOT_JS = Bonito.ES6Module(joinpath(@__DIR__, "TimePlots.js"))

Base.showable(::MIME"text/html", ::TimePlot) = true
Base.show(io::IO, ::MIME"text/plain", ::TimePlot) = print(io, "TimePlot (live Plotly viewer; display as HTML)")
Base.show(io::IO, mime::MIME"text/html", plot::TimePlot) =
    show(io, mime, Bonito.App(session -> Bonito.jsrender(session, plot)))

function Bonito.jsrender(session::Bonito.Session, plot::TimePlot)
    template = plot_template(plot.source)
    visibility = [get(trace, :visible, true) for trace in template.data]
    interval = get(get(template.layout, :xaxis, Dict()), :range, [])
    interval = length(interval) == 2 ? vec(interval) : template.full_range
    initial = window_data(plot.source, interval, 1000, visibility)
    overview = window_data(plot.source, template.full_range, 1000, fill(true, length(visibility)))
    request = Bonito.Observable(Dict{String,Any}())
    response = Bonito.Observable(PlotlyBase.JSON.json((; initial..., id=0)))
    Bonito.on(session, request) do query
        isempty(query) && return
        result = window_data(plot.source, query["range"], query["width"], query["visibility"])
        response[] = PlotlyBase.JSON.json((; result..., id=query["id"]))
        return nothing
    end
    width = get(template.layout, :width, nothing)
    height = get(template.layout, :height, nothing)
    style = "width:$(isnothing(width) ? "100%" : "$(width)px");height:$(isnothing(height) ? "100%" : "$(height)px");min-height:360px;"
    root = Bonito.DOM.div(; class="koma-time-plot", style)
    Bonito.onload(session, root, js"""async root => {
        const [Plotly, viewer] = await Promise.all([$(TIME_PLOTLY), $(TIME_PLOT_JS)]);
        viewer.mount(root, Plotly, $(request), $(response), $(PlotlyBase.JSON.json(overview)));
    }""")
    return Bonito.jsrender(session, root)
end

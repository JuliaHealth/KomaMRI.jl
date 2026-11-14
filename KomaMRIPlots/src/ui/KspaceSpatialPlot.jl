const KSPACE_DISPLAY_BUDGET = 256 * 256
const KSPACE_NAVIGATOR_HEIGHT = 120

struct KspaceSpatialSource{K,A,T,S,O}
    trajectory::K
    adc::A
    times::T
    trajectory_times::S
    overview::O
end

function kspace_spatial_plot(seq; width, height, darkmode, view_2d, adaptive=true, adaptive_navigator=false)
    sampled = discretize(seq; sampling_rule=KomaMRIBase.MomentSamplingRule(KomaMRIBase.MaxStepSizeRule(1, 5e-5)))
    rf_idx, rf_types = KomaMRIBase.get_RF_types(seq, sampled.t[1:end-1])
    trajectory, adc = get_kspace(sampled; rf_idx, rf_types)
    times = sampled.t[sampled.ADC]
    extent = isempty(adc) ? trajectory : adc
    bounds = [minimum(extent; dims=1); maximum(extent; dims=1)]
    all(iszero, bounds[2, :] - bounds[1, :]) && (bounds .+= [-1; 1])
    template = kspace_plot(adaptive ? bounds : trajectory, adaptive ? bounds : adc,
        adaptive ? zeros(2) : times, length(times);
        width, height, darkmode, view_2d)
    if adaptive
        for trace in template.data[1:2], axis in (:x, :y, :z)
            haskey(trace, axis) && (trace[axis] = Float64[])
        end
        template.data[2][:text] = Float64[]
        template.data[2][:marker][:color] = Float64[]
    end
    navigator = if adaptive_navigator
        plot = plot_seq(seq; adaptive=true, slider=false, height=KSPACE_NAVIGATOR_HEIGHT,
            range=[0.0, dur(seq) * 1e3], darkmode)
        layout = plot_template(plot.source).layout
        layout[:showlegend] = false
        layout[:updatemenus] = []
        layout[:yaxis][:fixedrange] = true
        plot
    else
        nothing
    end
    source = KspaceSpatialSource(trajectory, adc, times, sampled.t[2:end], navigator)
    return adaptive ? SpatialPlot(source, template, view_2d ? [1, 2] : [1, 2, 3]) : pluto_plot(template, source)
end

spatial_controls(source::KspaceSpatialSource) = isnothing(source.overview) ? nothing :
    (; time_max=last(source.trajectory_times) * 1e3)

function spatial_navigator(source::KspaceSpatialSource)
    isnothing(source.overview) && return nothing
    height = get(plot_template(source.overview.source).layout, :height, 280)
    return Bonito.DOM.div(source.overview;
        style="width:100%;height:$(height)px;max-height:45%;flex:none;min-height:0;")
end

struct KspacePlot{S,P}
    source::S
    plot::P
end

kspace_controls(plot::KspacePlot) = (; spatial_controls(plot.source)...,
    times=(; trajectory=plot.source.trajectory_times .* 1e3, adc=plot.source.times .* 1e3))

Base.show(io::IO, ::MIME"text/plain", ::KspacePlot) = print(io, "KspacePlot (live sequence navigator; display as HTML)")

function Bonito.jsrender(session::Bonito.Session, plot::KspacePlot)
    width, height = get(plot.plot.layout, :width, nothing), get(plot.plot.layout, :height, nothing)
    style = "width:$(isnothing(width) ? "100%" : "$(width)px");height:$(isnothing(height) ? "600px" : "$(height)px");min-height:$(isnothing(height) ? 300 : 0)px;"
    root = Bonito.DOM.div(TIME_PLOTLY, SPATIAL_PLOT_JS, spatial_navigator(plot.source);
        style, dataLinkedTimePlots=true)
    figure, controls = PlotlyBase.JSON.json(plot.plot), PlotlyBase.JSON.json(kspace_controls(plot))
    Bonito.onload(session, root, js"""async root => {
        const [Plotly, viewer] = await Promise.all([$(TIME_PLOTLY), $(SPATIAL_PLOT_JS)]);
        const graph = document.createElement('div');
        root.prepend(graph);
        const figure = JSON.parse($(figure));
        await Plotly.newPlot(graph, figure.data, figure.layout, figure.config);
        viewer.mountKspace(root, graph, Plotly, JSON.parse($(controls)));
    }""")
    return Bonito.jsrender(session, root)
end

function kspace_selection(source, query)
    adc = eachindex(source.times)
    interval = (0.0, last(source.trajectory_times))
    if haskey(query, "time_range")
        interval = Float64.(query["time_range"]) ./ 1e3
        adc = searchsortedfirst(source.times, first(interval)):searchsortedlast(source.times, last(interval))
    end
    vertices = searchsortedfirst(source.trajectory_times, first(interval)):searchsortedlast(source.trajectory_times, last(interval))
    return adc, first(vertices):last(vertices)-1
end

# Clip a projected segment, including crossings with both endpoints outside the view.
function visible_segment(a, b)
    all(isfinite, a) && all(isfinite, b) || return nothing
    lo, hi = 0.0, 1.0
    for d in 1:3
        delta = b[d] - a[d]
        if iszero(delta)
            -1 <= a[d] <= 1 || return nothing
        else
            first, last = minmax((-1-a[d]) / delta, (1-a[d]) / delta)
            lo, hi = max(lo, first), min(hi, last)
            lo <= hi || return nothing
        end
    end
    middle = (lo + hi) / 2
    return ntuple(d -> a[d] + middle * (b[d] - a[d]), 3)
end

function kspace_segments(trajectory, matrix, width, height, budget; selection=1:size(trajectory, 1)-1)
    isempty(selection) && return Int[], 0
    nx = clamp(floor(Int, sqrt(budget * width / height)), 1, budget)
    ny = max(1, budget ÷ nx)
    cells = Dict{Tuple{Int,Int},Int}()
    visible = Int[]
    available = 0
    previous = projected_point(matrix, @view trajectory[first(selection), :])
    for i in selection
        next = projected_point(matrix, @view trajectory[i+1, :])
        middle = visible_segment(previous, next)
        previous = next
        isnothing(middle) && continue
        available += 1
        available <= budget && push!(visible, i)
        cell = (clamp(floor(Int, (middle[1]+1) * nx / 2), 0, nx-1),
            clamp(floor(Int, (middle[2]+1) * ny / 2), 0, ny-1))
        old = get(cells, cell, i)
        cells[cell] = hash(i) < hash(old) ? i : old
    end
    return available <= budget ? visible : sort!(collect(values(cells))), available
end

function spatial_data(plot::SpatialPlot{<:KspaceSpatialSource}, query)
    source = plot.source
    matrix = reshape(Float64.(query["projection"]), 4, 4)
    width, height = max(1.0, Float64(query["width"])), max(1.0, Float64(query["height"]))
    budget = clamp(floor(Int, width * height / SPATIAL_PIXEL_SIZE^2), 1, KSPACE_DISPLAY_BUDGET)
    visible = get(query, "visible", [true, true, true])
    selected_adc, selected_segments = kspace_selection(source, query)
    segments, available_segments = visible[1] ?
        kspace_segments(source.trajectory, matrix, width, height, budget; selection=selected_segments) : (Int[], 0)
    adc_budget = length(plot.dimensions) == 2 ? max(1, budget ÷ 2) : budget
    scale = sqrt(adc_budget * SPATIAL_PIXEL_SIZE^2 / (width * height))
    coordinates = ntuple(d -> @view(source.adc[:, d]), 3)
    ids, available_adc = visible[2] ? spatial_indices(coordinates, eachindex(source.times), matrix,
        width * scale, height * scale, length(plot.dimensions); selection=selected_adc) : (Int[], 0)
    data = deepcopy(plot.template.data)
    # Join consecutive segments; only omitted segments need a line break.
    vertices = Int[]
    for i in segments
        if isempty(vertices) || last(vertices) != i
            isempty(vertices) || push!(vertices, 0)
            push!(vertices, i)
        end
        push!(vertices, i+1)
    end
    for (axis, d) in zip((:x, :y, :z), plot.dimensions)
        data[1][axis] = [iszero(i) ? nothing : source.trajectory[i, d] for i in vertices]
        data[2][axis] = source.adc[ids, d]
    end
    data[1][:customdata] = vertices
    data[2][:customdata] = ids
    data[2][:text] = round.(source.times[ids] .* 1e3; digits=3)
    data[2][:marker][:color] = ids .- 1
    shown = length(segments) + length(ids)
    available = available_segments + available_adc
    trajectory_detail = visible[1] ? "$(length(segments)) / $available_segments visible segments" : "hidden"
    adc_detail = visible[2] ? "$(length(ids)) / $available_adc visible samples" : "hidden"
    return (; data, shown, available, coarse=false,
        detail="Trajectory: $trajectory_detail; ADC: $adc_detail.")
end

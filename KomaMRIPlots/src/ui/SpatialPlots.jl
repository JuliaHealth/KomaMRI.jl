const SPATIAL_PIXEL_SIZE = 2
const SPATIAL_FRAME_CACHE = 2
const SPATIAL_MOTION_BUDGET = 10_000

"""
    plot_phantom(obj; key=:ρ, adaptive=false, kwargs...)

Plot intrinsic phantom properties and motion-derived speed/velocity in one viewer.
Property buttons share spin coordinates, camera and the time slider. `kwargs` are
the display options of [`plot_phantom_map`](@ref). Velocities are in cm/s; magenta
and cyan mark exit and re-entry, excluding reset jumps from velocity estimates.

The default is a self-contained Plotly plot. With `adaptive=true`, Julia supplies
visible spins on demand and playback loops over the interval in five seconds,
skipping display frames when necessary; every time step remains selectable.

# Keywords — both modes
- `key=:ρ`: initially selected property.
- `time_samples=0`: add an evenly spaced time grid to the native motion timepoints;
  duplicates are removed, and native timepoints are never discarded. Static phantoms
  have one frame. Nonadaptive plots embed every resulting frame.
- `view_2d`: defaults to true for phantoms with fewer than three spatial dimensions.
- `height`, `width`, `darkmode`, `colorbar`, and color limits: see [`plot_phantom_map`](@ref).

# Keywords — only `adaptive=false`
- `max_spins`: cap the embedded spin subset; defaults to 20,000 when omitted.
  Passing it with `adaptive=true` throws `ArgumentError`; live sampling uses the viewport.
"""
function plot_phantom(obj; key=:ρ, adaptive=false, max_spins=nothing,
    time_samples=0, view_2d=sum(KomaMRIBase.get_dims(obj)) < 3, kwargs...)
    return phantom_plot(obj, key; adaptive, max_spins, time_samples, view_2d,
        properties=nothing, kwargs...)
end

phantom_properties(::NoMotion) = (:ρ, :T1, :T2, :T2s, :Δw)
phantom_properties(motion) = (:ρ, :T1, :T2, :T2s, :Δw, :speed, :vx, :vy, :vz)

"""
    SpatialPlot

Live Plotly spatial viewer returned by `plot_phantom_map(...; adaptive=true)` and
`plot_coil_sens(...; adaptive=true)` and `plot_kspace(...; adaptive=true)`.
Zooming, rotating and changing time or coil
fetch detail from Julia; keep the Julia session running, including in notebooks.
The default `adaptive=false` remains a self-contained `PlotlyBase.Plot`.

Phantom sampling retains original spin IDs and property values. Two-dimensional
views retain screen-cell extrema; 3D views sample throughout spatial voxels.
All spins remain available on zoom. Only spins inside the viewport consume its
sampling budget, including during motion previews.
`max_spins` is only accepted with `adaptive=false`. Keyword availability is documented
on each `plot_*` function. Motion bounds are scanned once, one frame at a time,
to keep axes fixed. The coordinate cache retains only two frames.
Camera and color limits remain unchanged when selecting a time.
Dragging or playing motion uses a stable spin preview; stopping restores spatial detail.
The final frame shows the recorded endpoint; playback wrapping is not a flow-reset event.

Coil maps refine a grid down to `spacing`, evaluating `get_sens` only at selected
positions. Their magnitude scale uses a model-wide bound shared by all coils.
An in-plot indicator marks reduced spatial detail, not missing simulation data.
"""
struct SpatialPlot{S,P}
    source::S
    template::P
    dimensions::Vector{Int}
end

struct PhantomSpatialSource{O,V,T,C,P}
    obj::O
    values::V
    times::T
    frames::C
    properties::P
end

struct CoilSpatialSource{R,G}
    receiver::R
    grid::G
end

struct PhantomVelocity
    component::Symbol
end

flow_resets!(reset, ::NoMotion, interval) = nothing
function flow_resets!(reset, motions::MotionList, interval)
    foreach(m -> flow_resets!(reset, m, interval), motions.motions)
    return nothing
end
function flow_resets!(reset, motion::Motion, interval)
    selected = @view reset[KomaMRIBase.get_indexing_range(motion.spins)]
    flow_resets!(selected, motion.action, motion.time, interval)
    return nothing
end
flow_resets!(reset, action, time, interval) = nothing
function flow_resets!(reset, action::FlowPath, time, interval)
    time = phantom_endpoint(time, maximum(interval))
    unit_interval = KomaMRIBase.unit_time(collect(interval), time)
    lower, upper = extrema(unit_interval .* (size(action.spin_reset, 2) - 1))
    # A reset column ends a discontinuous interpolation interval, not a physical velocity.
    for column in floor(Int, lower)+2:ceil(Int, upper)+1
        reset .|= @view action.spin_reset[:, column]
    end
    return nothing
end

phantom_endpoint(motion::NoMotion, t) = motion
function phantom_endpoint(time::TimeCurve, t)
    (!time.periodic || t != last(times(time))) && return time
    return TimeCurve(; time.t, time.t_unit, time.periods, time.t_start, time.t_end, periodic=false)
end
phantom_endpoint(motion::Motion, t) =
    Motion(motion.action, phantom_endpoint(motion.time, t), motion.spins)
phantom_endpoint(motion::MotionList, t) = MotionList(map(m -> phantom_endpoint(m, t), motion.motions))

phantom_frame_values(values, source, frame, coordinates) = (values, nothing)
function phantom_frame_values(velocity::PhantomVelocity, source, frame, coordinates)
    length(source.times) == 1 && return (zeros(length(source.obj)), falses(length(source.obj)))
    previous = frame == 1 ? 2 : frame - 1
    other = phantom_frame(source, previous)
    interval = (source.times[previous], source.times[frame])
    dt = interval[2] - interval[1]
    values = if velocity.component == :speed
        sqrt.(sum((a .- b).^2 for (a, b) in zip(coordinates, other))) ./ abs(dt)
    else
        axis = findfirst(==(velocity.component), (:vx, :vy, :vz))
        (coordinates[axis] .- other[axis]) ./ dt
    end
    reset = falses(length(values))
    flow_resets!(reset, source.obj.motion, interval)
    values[reset] .= 0
    return values, reset
end

spatial_times(::NoMotion, count) = [0.0]
function spatial_times(motion, count)
    nodes = times(motion)
    return sort!(unique([nodes; range(extrema(nodes)...; length=max(count, 2))]))
end

phantom_coordinates(motion, obj, t) = get_spin_coords(motion, obj.x, obj.y, obj.z, [t]')
phantom_coordinates(motion::Motion, obj, t) = phantom_coordinates(motion.action, motion, obj, t)
phantom_coordinates(action, motion, obj, t) = get_spin_coords(motion, obj.x, obj.y, obj.z, [t]')
function phantom_coordinates(action::KomaMRIBase.ArbitraryAction, motion, obj, t)
    unit_t = only(KomaMRIBase.unit_time([t]', motion.time))
    knots = collect(range(zero(unit_t), one(unit_t); length=size(action.dx, 2)))
    column = clamp(searchsortedfirst(knots, unit_t) - 1, 1, length(knots) - 1)
    weight = (unit_t - knots[column]) / (knots[column + 1] - knots[column])
    selected = KomaMRIBase.get_indexing_range(motion.spins)
    # A frame queries every spin at the same time; only the time axis needs interpolation.
    return map((obj.x, obj.y, obj.z), (action.dx, action.dy, action.dz)) do position, displacement
        coordinates = position .+ zero(eltype(position)) .* [t]'
        @views coordinates[selected, :] .+= (1 - weight) .* displacement[:, column] .+
            weight .* displacement[:, column + 1]
        coordinates
    end
end

function phantom_frame(source, frame)
    hit = findfirst(pair -> first(pair) == frame, source.frames)
    !isnothing(hit) && return last(source.frames[hit])
    obj = source.obj
    t = source.times[frame]
    # The final recorded sample is the end of the trajectory, not the next playback cycle.
    motion = frame == length(source.times) ? phantom_endpoint(obj.motion, t) : obj.motion
    coordinates = map(v -> 100 .* vec(v), phantom_coordinates(motion, obj, t))
    length(source.frames) == SPATIAL_FRAME_CACHE && popfirst!(source.frames)
    push!(source.frames, frame => coordinates)
    return coordinates
end

phantom_units(key) = key in (:T1, :T2, :T2s) ? (1e3, " ms") :
    key in (:x, :y, :z) ? (100, " cm") : key == :Δw ? (1 / 2π, " Hz") : (1, "")

function phantom_property_plot(obj, key; view_2d, height=700, width=nothing,
    darkmode=false, colorbar=true, kwargs...)
    velocity = key in (:speed, :vx, :vy, :vz)
    label = key == :speed ? "|v|" : string(key)
    factor, unit = velocity ? (1, " cm/s") : phantom_units(key)
    cmin, cmax = velocity ? (0.0, 1.0) : extrema(getproperty(obj, key)) .* factor
    key ∉ (:x, :y, :z, :Δw) && (cmin = 0)
    key == :T1 && (cmax = 2500)
    key == :T2 && (cmax = 250)
    colorscale = if velocity
        key == :speed ? "Viridis" : [[0, "blue"], [0.5, "white"], [1, "red"]]
    elseif key in (:T1, :T2, :T2s)
        colors = replace.(string.(relaxationColorMap(key == :T1 ? "T1" : "T2") .* 255),
            "RGB{Float64}" => "rgb")
        collect(zip(range(0.0, 1.0; length=length(colors)), colors))
    else
        "Greys"
    end
    cmin, cmax = get(kwargs, :zmin, cmin), get(kwargs, :zmax, cmax)
    dimensions = findall(KomaMRIBase.get_dims(obj))
    isempty(dimensions) && (dimensions = [1, 2])
    length(dimensions) == 1 && (dimensions = sort([only(dimensions), mod1(only(dimensions) + 1, 3)]))
    background, text, plot_background, grid_color, _ = theme_chooser(darkmode)
    axes = [attr(; title=name, ticksuffix=" cm", backgroundcolor=plot_background,
        gridcolor=grid_color, zerolinecolor=grid_color) for name in ("x", "y", "z")]
    layout = Layout(; title=obj.name * ": " * string(key), paper_bgcolor=background,
        plot_bgcolor=plot_background, font_color=text,
        margin=attr(; t=50, l=0, r=0, b=view_2d ? 50 : 0),
        modebar=attr(; orientation="h", bgcolor=background, color=text, activecolor=plot_background))
    if view_2d
        layout[:xaxis], layout[:yaxis] = axes[dimensions[1]], axes[dimensions[2]]
        layout[:xaxis][:scaleanchor], layout[:yaxis][:scaleratio] = "y", 1
        layout[:xaxis][:constrain] = layout[:yaxis][:constrain] = "domain"
        layout[:autosize] = true
    else
        foreach(axis -> axis[:fixedrange] = true, axes)
        layout[:scene] = attr(; xaxis=axes[1], yaxis=axes[2], zaxis=axes[3],
            aspectmode="manual", aspectratio=attr(; x=1, y=1, z=1))
    end
    hover = "x: %{x:.1f} cm<br>y: %{y:.1f} cm<br>" * (view_2d ? "" : "z: %{z:.1f} cm<br>")
    trace = (view_2d ? scattergl : scatter3d)(; x=Float64[], y=Float64[],
        (view_2d ? (;) : (; z=Float64[]))..., mode="markers", visible=true, showlegend=false,
        marker=attr(; color=Float64[], showscale=colorbar, colorscale,
            colorbar=attr(; ticksuffix=unit, title=label), cmin, cmax, size=view_2d ? 4 : 2),
        text=Float64[], hovertemplate=hover * "<b>$label</b>: %{text}$unit<extra></extra>")
    height !== nothing && (layout.height = height)
    width !== nothing && (layout.width = width)
    config = PlotConfig(; displaylogo=false, toImageButtonOptions=attr(; format="svg").fields,
        modeBarButtonsToRemove=["zoom", "pan", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"])
    return PlotlyBase.Plot([trace], layout; config)
end

function phantom_velocity_limits!(limits, coordinates, previous, dt, reset)
    for spin in eachindex(reset)
        reset[spin] && continue
        vx, vy, vz = map((a, b) -> (a[spin] - b[spin]) / dt, coordinates, previous)
        maxima = (sqrt(vx^2 + vy^2 + vz^2), abs(vx), abs(vy), abs(vz))
        for i in eachindex(limits)
            limits[i] = max(limits[i], maxima[i])
        end
    end
    return nothing
end

function phantom_spatial_plot(obj, key; frame_times, view_2d, properties, kwargs...)
    property_plots = NamedTuple(property => phantom_property_plot(obj, property; view_2d, kwargs...)
        for property in properties)
    template = deepcopy(property_plots[key])
    first_frame = map(v -> 100 .* vec(v), phantom_coordinates(obj.motion, obj, first(frame_times)))
    values = NamedTuple(property => begin
        property in (:speed, :vx, :vy, :vz) ? PhantomVelocity(property) :
            getproperty(obj, property) .* first(phantom_units(property))
    end for property in properties)
    source = PhantomSpatialSource(obj, values, frame_times, [1 => first_frame], property_plots)
    extent = maximum(v -> maximum(abs, v), first_frame)
    bounds = map(extrema, first_frame)
    velocity_limits = zeros(4)
    has_velocity = any(in((:speed, :vx, :vy, :vz)), properties)
    previous = first_frame
    reset = falses(length(obj))
    for frame in eachindex(frame_times)
        coordinates = phantom_frame(source, frame)
        extent = max(extent, maximum(v -> maximum(abs, v), coordinates))
        bounds = map((bound, v) -> (min(bound[1], minimum(v)), max(bound[2], maximum(v))), bounds, coordinates)
        if has_velocity && frame > 1
            interval = (frame_times[frame - 1], frame_times[frame])
            dt = interval[2] - interval[1]
            fill!(reset, false)
            flow_resets!(reset, obj.motion, interval)
            phantom_velocity_limits!(velocity_limits, coordinates, previous, dt, reset)
        end
        previous = coordinates
    end
    extent = iszero(extent) ? 1.0 : extent
    layout = template.layout
    for (i, property) in enumerate((:speed, :vx, :vy, :vz))
        property in properties || continue
        trace = first(property_plots[property].data)
        limit = get(kwargs, :zmax, iszero(velocity_limits[i]) ? 1.0 : velocity_limits[i])
        trace[:marker][:cmin] = property == :speed ? 0 : -limit
        trace[:marker][:cmax] = limit
    end
    template.data[1] = deepcopy(first(property_plots[key].data))
    if has_velocity
        layout[:legend] = attr(; orientation="h", xref="container", yref="container", yanchor="top", x=0.5, xanchor="center", y=1)
    end
    layout.fields[:title] = length(properties) > 1 ?
        attr(; text=obj.name, x=0.5, xanchor="center", y=1, yref="paper",
            yanchor="bottom", pad=attr(; b=16)) : obj.name * ": " * string(key)
    layout.fields[:meta] = Dict("phantom_property" => string(key), "spatial_range" => [-extent, extent], "spatial_bounds" => bounds)
    if length(properties) > 1
        groups = has_velocity ? (properties[1:end-4], properties[end-3:end]) : (properties,)
        layout[:updatemenus] = [attr(; name="phantom-property", type="buttons", direction="right",
            x=0.5, xanchor=has_velocity ? (group_index == 1 ? "right" : "left") : "center",
            y=1, yanchor="bottom", pad=attr(; b=has_velocity ? 95 : 50, l=has_velocity ? 10 : 0, r=has_velocity ? 10 : 0), showactive=true,
            bgcolor="#2a7fb8", bordercolor="#2a7fb8", font=attr(; color="white"),
            active=something(findfirst(==(key), group), 0)-1,
            buttons=[attr(; label=property == :speed ? "|v|" : string(property), method="skip", name="phantom-property",
                args=[string(property)]) for property in group]) for (group_index, group) in enumerate(groups)]
        has_velocity && (layout[:annotations] = [attr(; text="—", x=0.5, y=1,
            xref="paper", yref="paper", xanchor="center", yanchor="middle", yshift=111, showarrow=false)])
        layout[:margin][:t] = has_velocity ? 135 : 84
    end
    length(frame_times) > 1 && (layout[:margin][:b] = 50)
    dimensions = view_2d ? [findfirst(==(layout[axis][:title]), ("x", "y", "z"))
        for axis in (:xaxis, :yaxis)] : [1, 2, 3]
    for axis in (view_2d ? (layout[:xaxis], layout[:yaxis]) :
        (layout[:scene][:xaxis], layout[:scene][:yaxis], layout[:scene][:zaxis]))
        axis[:range] = [-extent, extent]
        axis[:autorange] = false
    end
    layout[:sliders] = [attr(; name="time", active=0, visible=length(frame_times)>1,
        pad=attr(; l=30, b=30), currentvalue_prefix="t = ", currentvalue_suffix=" ms",
        steps=[attr(; label=round(t * 1e3; digits=2), value=string(i), method="skip")
            for (i, t) in enumerate(frame_times)])]
    return SpatialPlot(source, template, dimensions)
end

function portable_phantom_frames(plot, query)
    data = GenericTrace[]
    for frame in eachindex(plot.source.times)
        traces = spatial_data(plot, merge(query, Dict("frame" => frame))).data
        for trace in traces
            trace[:visible] = frame == 1
            trace[:uid] *= "-$(bitstring(plot.source.times[frame]))"
        end
        append!(data, traces)
    end
    return data
end

function portable_phantom_plot(plot)
    source, template = plot.source, deepcopy(plot.template)
    query = spatial_query(plot; width=2length(source.obj), height=2)
    # All spins belong to the portable plot; camera clipping is handled by Plotly.
    query["projection"] = zeros(16)
    query["projection"][16] = 1
    initial_key = query["property"]
    data = portable_phantom_frames(plot, query)
    traces_per_frame = length(data) ÷ length(source.times)
    menus = get(template.layout, :updatemenus, [])
    for (menu_index, menu) in enumerate(menus), (button_index, button) in enumerate(menu[:buttons])
        query["property"] = only(button[:args])
        attributes = Dict(name => Any[] for name in ("text", "hovertemplate", "showlegend",
            "marker.color", "marker.colorscale", "marker.cmin", "marker.cmax",
            "marker.colorbar", "marker.size"))
        for frame in eachindex(source.times)
            query["frame"] = frame
            traces = query["property"] == initial_key ?
                @view(data[(frame-1)*traces_per_frame+1:frame*traces_per_frame]) : spatial_data(plot, query).data
            for trace in traces, (name, values) in attributes
                path = split(name, '.')
                push!(values, length(path) == 1 ? trace[Symbol(name)] : trace[:marker][Symbol(last(path))])
            end
        end
        active = Dict("updatemenus[$(i-1)].active" => i == menu_index ? button_index-1 : -1
            for i in eachindex(menus))
        button[:method], button[:args] = "update", [attributes, active]
    end
    for (frame, step) in enumerate(only(template.layout[:sliders])[:steps])
        step[:method] = "restyle"
        step[:args] = [attr(; visible=[cld(i, traces_per_frame) == frame for i in eachindex(data)])]
    end
    return PlotlyBase.Plot(data, template.layout; config=template.config)
end

coil_color_limit(::UniformCoilSens) = 1.0
coil_color_limit(receiver::ArbitraryCoilSens) = maximum(abs, receiver.coil_sens)
function coil_color_limit(receiver::BirdcageCoilSens)
    distance = (1 - BIRDCAGE_PLOT_RADIUS_FRACTION) * receiver.radius
    return 2receiver.L / (distance * hypot(distance, receiver.L))
end

function spatial_projection(plot)
    matrix = zeros(4, 4)
    matrix[4, 4] = 1
    if length(plot.dimensions) == 2
        for (row, axis) in enumerate((:xaxis, :yaxis))
            lower, upper = plot.template.layout[axis][:range]
            matrix[row, plot.dimensions[row]] = 2 / (upper - lower)
            matrix[row, 4] = -(upper + lower) / (upper - lower)
        end
    else
        lower, upper = plot.template.layout[:scene][:xaxis][:range]
        extent = (upper - lower) / 2
        matrix[1, 1:3] = [1, -1, 0] ./ (2extent)
        matrix[2, 1:3] = [1, 1, -2] ./ (3extent)
    end
    return matrix
end

function projected_point(matrix, position)
    p = ntuple(row -> sum(matrix[row, col] * position[col] for col in 1:3) + matrix[row, 4], 4)
    p[4] > 0 || return (Inf, Inf, Inf)
    return ntuple(i -> p[i] / p[4], 3)
end

function voxel_indices(coordinates, visible, budget)
    bounds = map(v -> extrema(@view v[visible]), coordinates)
    divisions = max(1, floor(Int, cbrt(budget)))
    scale = map(bounds) do (lower, upper)
        iszero(upper - lower) ? 0.0 : divisions / (upper - lower)
    end
    cells = Dict{NTuple{3,Int},Int}()
    for i in visible
        position = ntuple(d -> (coordinates[d][i] - bounds[d][1]) * scale[d], 3)
        cell = map(v -> min(floor(Int, v), divisions - 1), position)
        # Stable spin-ID selection avoids both intensity bias and a regular voxel-center lattice.
        previous = get(cells, cell, i)
        cells[cell] = hash(i) < hash(previous) ? i : previous
    end
    return sort!(collect(values(cells)))
end

function spatial_indices(coordinates, values, matrix, width, height, dimensions;
    preview=false, selection=eachindex(values))
    width, height = Float64(width), Float64(height)
    budget = max(1, floor(Int, width * height / SPATIAL_PIXEL_SIZE^2))
    visible = Int[]
    cells = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for i in selection
        x, y, z = projected_point(matrix, getindex.(coordinates, i))
        all(isfinite, (x, y, z)) && -1 <= x <= 1 && -1 <= y <= 1 && -1 <= z <= 1 || continue
        push!(visible, i)
        (dimensions == 3 || preview) && continue
        cell = (floor(Int, (x + 1) * width / (2SPATIAL_PIXEL_SIZE)),
            floor(Int, (y + 1) * height / (2SPATIAL_PIXEL_SIZE)))
        lo, hi = get(cells, cell, (i, i))
        cells[cell] = (values[i] < values[lo] ? i : lo, values[i] > values[hi] ? i : hi)
    end
    if preview
        indices = partialsort!(copy(visible), 1:min(length(visible), budget, SPATIAL_MOTION_BUDGET); by=hash)
        return sort!(indices), length(visible)
    end
    length(visible) <= budget && return visible, length(visible)
    dimensions == 3 && return voxel_indices(coordinates, visible, budget), length(visible)
    indices = sort!(unique!([i for pair in Base.values(cells) for i in pair]))
    return indices, length(visible)
end

function spatial_data(plot::SpatialPlot{<:PhantomSpatialSource}, query)
    haskey(plot.template.layout[:meta], "phantom_slices") && return phantom_slice_data(plot, query)
    return phantom_spatial_data(plot, query, plot.dimensions)
end

function phantom_spatial_data(plot, query, dimensions; selection=eachindex(plot.source.obj.x))
    preview = get(query, "preview", false)
    source = plot.source
    key = Symbol(get(query, "property", first(keys(source.values))))
    frame = Int(query["frame"])
    coordinates = phantom_frame(source, frame)
    values, reset = phantom_frame_values(source.values[key], source, frame, coordinates)
    # Browser arrays arrive as Vector{Any}; normalize before the per-spin work.
    matrix = reshape(Float64.(query["projection"]), 4, 4)
    selection_values = length(source.values) == 1 ? values : source.obj.ρ
    indices, available = spatial_indices(coordinates, selection_values, matrix,
        query["width"], query["height"], length(dimensions); preview, selection)
    shown = length(indices)
    velocity = key in (:speed, :vx, :vy, :vz)
    if isnothing(reset) && length(source.values) > 1
        reset = falses(length(values))
        frame > 1 && flow_resets!(reset, source.obj.motion, (source.times[frame - 1], source.times[frame]))
    end
    groups = if isnothing(reset)
        [indices]
    else
        exiting = falses(length(values))
        frame < length(source.times) && flow_resets!(exiting, source.obj.motion,
            (source.times[frame], source.times[frame + 1]))
        # The first frame uses a forward velocity difference, but has no preceding entry.
        frame == 1 && fill!(reset, false)
        [filter(i -> !reset[i] && !exiting[i], indices),
            filter(i -> exiting[i] && !reset[i], indices), filter(i -> reset[i], indices)]
    end
    data = map(enumerate(groups)) do (group, ids)
        trace = deepcopy(first(source.properties[key].data))
        trace[:uid] = "phantom-$group"
        for (axis, dim) in zip((:x, :y, :z), dimensions)
            trace[axis] = coordinates[dim][ids]
        end
        trace[:ids] = string.(ids)
        trace[:customdata] = ids
        trace[:marker][:color] = values[ids]
        trace[:text] = round.(values[ids]; digits=4)
        if group > 1
            name = group == 2 ? "Exit" : "Re-entry"
            trace[:name] = name
            trace[:showlegend] = velocity
            velocity && (trace[:marker][:color] = group == 2 ? "#ff40df" : "#00e5ff")
            velocity && (trace[:marker][:size] *= 2)
            trace[:marker][:showscale] = false
            velocity && (trace[:hovertemplate] = "Spin %{customdata}<br>$name<extra></extra>")
            # Empty event traces otherwise disappear from Plotly's legend between frames.
            if isempty(ids)
                for axis in (:x, :y, :z)[1:length(dimensions)]
                    trace[axis] = [nothing]
                end
            end
        end
        trace
    end
    return (; data, shown, available, coarse=false, frame, preview, property=string(key))
end

function phantom_plot(obj, key; adaptive, max_spins, time_samples,
    properties=(key,), kwargs...)
    adaptive && !isnothing(max_spins) && throw(ArgumentError("max_spins requires adaptive=false."))
    max_spins = something(max_spins, 20_000)
    if !adaptive && length(obj) > max_spins
        obj = obj[1:cld(length(obj), max_spins):end]
        @warn "For performance reasons, the number of displayed spins was capped to `max_spins`=$max_spins."
    end
    properties = isnothing(properties) ? phantom_properties(obj.motion) : properties
    frame_times = spatial_times(obj.motion, time_samples)
    plot = phantom_spatial_plot(obj, key; frame_times, properties, kwargs...)
    adaptive && return plot
    return pluto_plot(portable_phantom_plot(plot))
end

function grid_priority(grid, node, matrix, width, height)
    corners = [projected_point(matrix, ntuple(d -> 100grid[d][i[d]], 3))
        for i in Iterators.product(((first(r), last(r)) for r in node)...)]
    all(p -> !isfinite(first(p)), corners) && return -1.0
    any(p -> !all(isfinite, p), corners) && return max(width, height)
    bounds = ntuple(d -> extrema(p[d] for p in corners), 3)
    any(bound -> last(bound) < -1 || first(bound) > 1, bounds) && return -1.0
    return max((bounds[1][2] - bounds[1][1]) * width / 2,
        (bounds[2][2] - bounds[2][1]) * height / 2)
end

function spatial_grid(source, matrix, width, height)
    width, height = Float64(width), Float64(height)
    grid = source.grid
    budget = max(1, floor(Int, width * height / SPATIAL_PIXEL_SIZE^2))
    if prod(length, grid) <= budget
        positions = [p for p in Iterators.product(grid...) if inside_coil_plot(source.receiver, p)]
        return ntuple(d -> [p[d] for p in positions], 3), false
    end
    root = map(g -> 1:length(g), grid)
    priority(node) = grid_priority(grid, node, matrix, width, height)
    nodes = [(; range=root, priority=priority(root))]
    while length(nodes) < budget
        sort!(nodes; by=node -> -node.priority)
        refined = empty(nodes)
        splits = 0
        for node in nodes
            node.priority < 0 && continue
            if length(nodes) + splits < budget && node.priority > SPATIAL_PIXEL_SIZE &&
                any(r -> length(r) > 1, node.range)
                dim = argmax(map(r -> length(r), node.range))
                lower, upper = first(node.range[dim]), last(node.range[dim])
                middle = (lower + upper) ÷ 2
                for half in (lower:middle, middle+1:upper)
                    child = Base.setindex(node.range, half, dim)
                    score = priority(child)
                    score >= 0 && push!(refined, (; range=child, priority=score))
                end
                splits += 1
            else
                push!(refined, node)
            end
        end
        nodes = refined
        splits == 0 && break
    end
    positions = [ntuple(d -> grid[d][(first(node.range[d]) + last(node.range[d])) ÷ 2], 3)
        for node in nodes]
    filter!(p -> inside_coil_plot(source.receiver, p), positions)
    coordinates = ntuple(d -> [p[d] for p in positions], 3)
    return coordinates, any(node -> prod(length, node.range) > 1, nodes)
end

function spatial_data(plot::SpatialPlot{<:CoilSpatialSource}, query)
    source = plot.source
    matrix = reshape(Float64.(query["projection"]), 4, 4)
    coordinates, coarse = spatial_grid(source, matrix, query["width"], query["height"])
    values = get_sens(source.receiver, coordinates...)[:, query["coil"]]
    phase = query["component"] == 2
    trace = deepcopy(first(plot.template.data))
    trace[:uid] = "coil-map"
    for (axis, dim) in zip((:x, :y, :z), plot.dimensions)
        trace[axis] = 100 .* coordinates[dim]
    end
    trace[:name] = "Coil $(query["coil"])"
    trace[:marker][:color] = phase ? angle.(values) : abs.(values)
    phase && (trace[:hovertemplate] = replace(trace[:hovertemplate], "|S|:" => "Phase (rad):"))
    coloraxis = deepcopy(plot.template.layout[:coloraxis])
    if phase
        coloraxis[:cmin], coloraxis[:cmax], coloraxis[:colorscale] = -π, π, PHASE_COLORSCALE
        coloraxis[:colorbar] = attr(; title="rad", thickness=16,
            tickvals=[-π, -π/2, 0, π/2, π], ticktext=["−π", "−π/2", "0", "π/2", "π"])
    end
    return (; data=[trace], coloraxis, coarse, shown=length(values), available=length(values),
        spacing=step(first(source.grid)))
end

function spatial_query(plot; width=900, height=600)
    query = Dict{String,Any}(
    "id"=>0, "projection"=>vec(spatial_projection(plot)), "width"=>width, "height"=>height,
    "frame"=>1, "coil"=>1, "component"=>1,
    "property"=>get(get(plot.template.layout, :meta, Dict()), "phantom_property", nothing))
    slices = get(get(plot.template.layout, :meta, Dict()), "phantom_slices", nothing)
    if !isnothing(slices)
        query["slices"] = Dict("center" => slices.center, "thickness" => slices.thickness,
            "ranges" => [[plot.template.layout[Symbol(axis)][:range] for axis in plane.axes] for plane in slices.planes],
            "sizes" => [[width * (plane.domain[1][2] - plane.domain[1][1]),
                height * (plane.domain[2][2] - plane.domain[2][1])] for plane in slices.planes])
    end
    return query
end

const SPATIAL_PLOT_JS = Bonito.ES6Module(joinpath(@__DIR__, "SpatialPlots.js"))
spatial_controls(source) = nothing
spatial_navigator(source) = nothing
Base.show(io::IO, ::MIME"text/plain", ::SpatialPlot) = print(io, "SpatialPlot (live Plotly viewer; display as HTML)")

function plot_payload(plot::SpatialPlot)
    initial = spatial_data(plot, spatial_query(plot))
    query = q -> PlotlyBase.JSON.json((; spatial_data(plot, q)..., id=q["id"]))
    template = PlotlyBase.JSON.json((; layout=plot.template.layout, config=plot.template.config,
        plot.dimensions, kspace=spatial_controls(plot.source)))
    return (; initial=PlotlyBase.JSON.json((; initial..., id=0)), template, query,
        layout=plot.template.layout)
end

function Bonito.jsrender(session::Bonito.Session, plot::SpatialPlot)
    payload = plot_payload(plot)
    request, response = plot_observables(session, payload)
    width, height = get(plot.template.layout, :width, nothing), get(plot.template.layout, :height, nothing)
    style = "width:$(isnothing(width) ? "100%" : "$(width)px");height:$(isnothing(height) ? "100%" : "$(height)px");min-height:$(isnothing(height) ? 300 : 0)px;"
    root = Bonito.DOM.div(TIME_PLOTLY, SPATIAL_PLOT_JS, spatial_navigator(plot.source);
        class="koma-spatial-plot", style, dataLinkedTimePlots=true)
    Bonito.onload(session, root, js"""async root => {
        const [Plotly, viewer] = await Promise.all([$(TIME_PLOTLY), $(SPATIAL_PLOT_JS)]);
        if (!root?.isConnected) return;
        viewer.mount(root, Plotly, $(request), $(response), $(Bonito.JSString(PlotlyBase.JSON.json(payload.template))));
    }""")
    return Bonito.jsrender(session, root)
end

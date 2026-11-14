const PHANTOM_SLICE_STEPS = 101
const PHANTOM_SLICE_PLANES = (
    (; name="Axial (x'-y)", dimensions=(1, 2), normal=3, axes=("xaxis", "yaxis"), refs=("x", "y"), domain=([0, 0.48], [0.535, 1])),
    (; name="Sagittal (z-y)", dimensions=(3, 2), normal=1, axes=("xaxis2", "yaxis2"), refs=("x2", "y2"), domain=([0.52, 1], [0.535, 1])),
    (; name="Coronal (x'-z)", dimensions=(1, 3), normal=2, axes=("xaxis3", "yaxis3"), refs=("x3", "y3"), domain=([0, 0.48], [0, 0.465])),
)
const SLICE_PLOT_JS = Bonito.ES6Module(joinpath(@__DIR__, "PhantomSlices.js"))

function slice_layout!(label, layout, planes, positions, center, ranges, spacing; reversed=(), scale=1)
    base = deepcopy(layout[:xaxis])
    for plane in planes, (row, name) in enumerate(plane.axes)
        dimension = plane.dimensions[row]
        axis = deepcopy(base)
        axis[:title], axis[:showticklabels] = "", false
        axis[:domain], axis[:anchor] = plane.domain[row], plane.refs[3-row]
        axis[:range] = dimension in reversed ? reverse(ranges[dimension]) : ranges[dimension]
        axis[:autorange], axis[:constrain] = false, "domain"
        delete!(axis, :scaleanchor)
        if row == 2
            axis[:scaleanchor] = plane.refs[1]
            axis[:scaleratio] = spacing[dimension] / spacing[plane.dimensions[1]]
        end
        layout[Symbol(name)] = axis
    end
    layout[:annotations] = [get(layout, :annotations, []); [attr(; name="slice-title-$(plane.normal)",
        text=plane.name, xref="$(plane.refs[1]) domain", yref="$(plane.refs[2]) domain",
        x=0, xanchor="left", y=1, yshift=18, showarrow=false, borderpad=0) for plane in planes];
        [attr(; name="slice-position-$(plane.normal)", text=label(plane.normal, center[plane.normal]),
        xref="$(plane.refs[1]) domain", yref="$(plane.refs[2]) domain",
        x=1, xanchor="right", y=1, yshift=18, showarrow=false, borderpad=0) for plane in planes]]
    layout[:sliders] = [attr(; name="slice-$(plane.normal)",
        active=searchsortedfirst(positions[plane.normal], center[plane.normal])-1,
        x=plane.domain[1][1], len=plane.domain[1][2]-plane.domain[1][1],
        y=plane.domain[2][1], yanchor="top", pad=attr(; t=6,l=0,r=0,b=0), font=attr(; size=10),
        ticklen=0, minorticklen=0, currentvalue=attr(; visible=false),
        steps=[attr(; label=scale==1 ? v : round(scale*v;digits=2),value=string(v),method="skip") for v in positions[plane.normal]])
        for plane in planes]
    return nothing
end

"""
    plot_phantom_slice(obj; thickness=5e-3, center=nothing, key=:ρ, adaptive=false, kwargs...)

Show linked axial (x–y), sagittal (z–y), and coronal (x–z) spin projections.
Coronal and axial reverse x.
`thickness` and the three-component `center` are in metres. The default center
uses spins nearest the middle of each coordinate range. Scroll a panel or use its
position slider to move its slab; click a spin to position the linked edge ticks.
Slabs are fixed in space and use current spin positions during motion.

Property and time controls are shared. `adaptive=true` filters slabs before
viewport sampling and requires live Julia. The default returns a self-contained
`PlotlyBase.Plot`, embedding all spins and motion frames so every slab remains
available offline. `max_spins` optionally caps that embedded set; it is not
accepted with `adaptive=true`. Other display keywords match [`plot_phantom`](@ref).
"""
function plot_phantom_slice(obj; thickness=5e-3, center=nothing, key=:ρ,
    adaptive=false, max_spins=nothing, time_samples=0, height=800, kwargs...)
    isfinite(thickness) && thickness > 0 || throw(ArgumentError("thickness must be finite and positive."))
    adaptive && !isnothing(max_spins) && throw(ArgumentError("max_spins requires adaptive=false."))
    if !isnothing(max_spins) && length(obj) > max_spins
        max_spins > 0 || throw(ArgumentError("max_spins must be positive."))
        obj = obj[1:cld(length(obj), max_spins):end]
        @warn "Displayed spins were capped to max_spins=$max_spins."
    end
    plot = phantom_spatial_plot(obj, key; frame_times=spatial_times(obj.motion, time_samples),
        view_2d=true, properties=phantom_properties(obj.motion), height, kwargs...)
    append!(plot.template.config.modeBarButtonsToRemove, ["select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d"])
    layout = plot.template.layout
    layout[:title] = ""
    layout[:showlegend] = false
    coordinates = phantom_frame(plot.source, 1)
    position = isnothing(center) ? map(coordinates) do v
        middle = (minimum(v) + maximum(v)) / 2
        v[argmin(abs.(v .- middle))]
    end : 100 .* collect(center)
    length(position) == 3 && all(isfinite, position) || throw(ArgumentError("center must contain three finite coordinates."))
    bounds = layout[:meta]["spatial_bounds"]
    model = (; planes=PHANTOM_SLICE_PLANES, center=position, thickness=100thickness, bounds)
    layout[:meta]["phantom_slices"] = model
    positions = map(1:3) do axis
        values = sort!(unique(coordinates[axis]))
        indices = round.(Int, range(1, length(values); length=min(length(values), PHANTOM_SLICE_STEPS)))
        sort!(unique!([values[indices]; bounds[axis]...; position[axis]]))
    end
    extent = layout[:meta]["spatial_range"]
    padding = 0.05 * (last(extent) - first(extent))
    extent = extent + [-padding, padding]
    layout[:meta]["spatial_range"] = extent
    control_padding = minimum(menu[:pad][:b] for menu in layout[:updatemenus]) - 24
    layout[:margin][:t] -= control_padding
    for menu in layout[:updatemenus]
        menu[:pad][:b] -= control_padding
    end
    time_slider = only(layout[:sliders])
    time_slider[:pad][:t] = 26
    time_slider[:pad][:l] = 55
    time_slider[:len] = 1
    time_slider[:currentvalue][:font] = attr(; size=get(layout[:font], :size, 12))
    time_slider[:font] = attr(; size=10)
    slice_layout!(layout, PHANTOM_SLICE_PLANES, positions, position, fill(extent,3), ones(3); reversed=(1,),scale=10) do axis,value
        "$(('x','y','z')[axis]) = $(round(10value; digits=2)) mm"
    end
    push!(layout[:sliders], time_slider)
    push!(layout[:updatemenus], attr(; name="slab-thickness", type="dropdown", direction="down",
        x=0.5, xanchor="center", y=1, yanchor="bottom",
        pad=attr(; b=minimum(menu[:pad][:b] for menu in layout[:updatemenus])), active=2, showactive=false,
        bgcolor="#2a7fb8", bordercolor="#2a7fb8", font=attr(; color="white"),
        buttons=[attr(; name="slab-thickness", label="$(round(1000thickness * factor; digits=2)) mm",
            method="skip", args=[100thickness * factor]) for factor in (0.25, 0.5, 1, 2, 4)]))
    layout[:margin][:l] = 20
    layout[:margin][:r] = 20
    layout.fields[:margin] = Dict{Symbol,Any}(layout[:margin])
    layout[:margin][:autoexpand] = false
    layout[:margin][:b] = length(plot.source.times) > 1 ? 104 : 38
    for property in values(plot.source.properties)
        first(property.data)[:marker][:colorbar][:len] = 0.465
        first(property.data)[:marker][:colorbar][:x] = 0.76
        first(property.data)[:marker][:colorbar][:xanchor] = "center"
        first(property.data)[:marker][:colorbar][:y] = 0.2325
    end
    initial = phantom_slice_data(plot, spatial_query(plot))
    plot.template.data = initial.data
    adaptive && return plot
    return portable_phantom_slice(plot)
end

function phantom_slice_data(plot, query)
    model = plot.template.layout[:meta]["phantom_slices"]
    slices = query["slices"]
    coordinates = phantom_frame(plot.source, Int(query["frame"]))
    panels = map(enumerate(model.planes)) do (panel, plane)
        center, thickness = Float64.(slices["center"]), Float64(slices["thickness"])
        selected = findall(v -> abs(v - center[plane.normal]) <= thickness / 2, coordinates[plane.normal])
        matrix = zeros(4, 4)
        matrix[4, 4] = 1
        for row in 1:2
            lo, hi = slices["ranges"][panel][row]
            matrix[row, plane.dimensions[row]] = 2 / (hi - lo)
            matrix[row, 4] = -(hi + lo) / (hi - lo)
        end
        width, height = slices["sizes"][panel]
        data = phantom_spatial_data(plot, merge(query, Dict("projection"=>vec(matrix), "width"=>width, "height"=>height)),
            plane.dimensions; selection=selected)
        for (group, trace) in enumerate(data.data)
            trace[:uid] = "phantom-slice-$panel-$group"
            trace[:xaxis], trace[:yaxis] = plane.refs
            trace[:marker][:showscale] = false
            trace[:showlegend] &= panel == 3
            colorbar = trace[:marker][:colorbar]
            trace[:hovertemplate] = "Spin %{customdata}<br>$(('x','y','z')[plane.dimensions[1]]): %{x:.2f} cm<br>" *
                "$(('x','y','z')[plane.dimensions[2]]): %{y:.2f} cm<br>" *
                "<b>$(colorbar[:title])</b>: %{text}$(colorbar[:ticksuffix])<extra></extra>"
        end
        data
    end
    # The scale must survive even when every slab is empty.
    marker = deepcopy(first(plot.source.properties[Symbol(query["property"])].data)[:marker])
    marker[:colorbar][:title] = attr(; text=marker[:colorbar][:title])
    marker[:color] = [marker[:cmin], marker[:cmax]]
    marker[:size], marker[:opacity] = 0, 0
    colorbar = scatter(; x=[0, 0], y=[0, 0], mode="markers", marker,
        uid="phantom-slice-colorbar", showlegend=false, hoverinfo="skip")
    return (; panels[1]..., data=vcat((panel.data for panel in panels)..., [colorbar]),
        shown=sum(panel.shown for panel in panels), available=sum(panel.available for panel in panels),
        slice_counts=[(; panel.shown, panel.available) for panel in panels])
end

function portable_phantom_slice(plot)
    source = plot.source
    layout = deepcopy(plot.template.layout)
    frames = map(eachindex(source.times)) do frame
        coordinates = phantom_frame(source, frame)
        entry, exit = falses(length(source.obj)), falses(length(source.obj))
        frame > 1 && flow_resets!(entry, source.obj.motion, (source.times[frame-1], source.times[frame]))
        frame < length(source.times) && flow_resets!(exit, source.obj.motion, (source.times[frame], source.times[frame+1]))
        (; coordinates, entry=findall(entry), exit=findall(exit))
    end
    values = Dict(string(key)=>value for (key, value) in pairs(source.values) if !(key in (:speed, :vx, :vy, :vz)))
    properties = Dict(string(key)=>PlotlyBase.JSON.lower(first(value.data)) for (key, value) in pairs(source.properties))
    layout[:meta]["phantom_slices"] = (; layout[:meta]["phantom_slices"]..., frames, values, properties, times=source.times)
    return pluto_plot(PlotlyBase.Plot(deepcopy(plot.template.data), layout; config=plot.template.config))
end

const UNIFORM_COIL_PLOT_FOV = 0.256
const BIRDCAGE_PLOT_RADIUS_FRACTION = 0.9

coil_plot_bounds(::UniformCoilSens) =
    ntuple(_ -> (-UNIFORM_COIL_PLOT_FOV / 2, UNIFORM_COIL_PLOT_FOV / 2), 3)
coil_plot_bounds(receiver::BirdcageCoilSens) = (
    (-receiver.radius, receiver.radius),
    (-receiver.radius, receiver.radius),
    (-receiver.L, receiver.L),
)
coil_plot_bounds(receiver::ArbitraryCoilSens) = map(extrema, (receiver.x, receiver.y, receiver.z))

"""
    get_coil_sens_fov(receiver)

Return the default plotting extent `(Lx, Ly, Lz)` in metres for a receive model.
Uniform sensitivities use a 256 mm cube, Birdcage uses its diameter and full length,
and arbitrary sensitivities use their sampled coordinate spans. A zero span denotes
a plane. This is a display extent, not a physical sensitivity cutoff.
"""
get_coil_sens_fov(receiver) = map(bounds -> last(bounds) - first(bounds), coil_plot_bounds(receiver))

inside_coil_plot(::AbstractRFReceiveSystem, position) = true
inside_coil_plot(receiver::BirdcageCoilSens, position) =
    hypot(position[1], position[2]) < BIRDCAGE_PLOT_RADIUS_FRACTION * receiver.radius

"""
    plot_coil_sens(sys; fov=nothing, spacing=0.01, height=700, width=nothing,
        darkmode=false, adaptive=false)

Plot receive-coil sensitivities without a phantom. Coordinates and `spacing` are in
metres; axes are displayed in centimetres. The sampling grid is centred in the receiver's
bounds and spaced by `spacing`, independently of its field of view.

`fov=(Lx, Ly, Lz)` overrides the displayed extent in metres, centred on the model's
bounds; a zero extent selects a plane. When omitted, [`get_coil_sens_fov`](@ref)
provides the extent. The override does not change the receiver or its sensitivities;
arbitrary maps return zero outside their sampled support.

Uniform sensitivities use a 256 mm box. Birdcage sensitivities use the cylinder's
half-length `L` and 90% of its radius, avoiding the ideal-wire singularities. Sampled
sensitivities use their stored coordinate bounds. Planar maps are plotted in their plane;
volumetric maps use a 3D point cloud with equal spatial scales.

The self-contained plot has a coil selector above it (only for multiple coils) and
Magnitude/Phase buttons below. Magnitude uses a grayscale range shared by all coils;
phase uses a cyclic scale over `[-π, π]`. Values come directly from `get_sens`, without
spin-density weighting or per-coil normalization.

With `adaptive=true`, return a live [`SpatialPlot`](@ref): refine the visible grid
on zoom or camera changes, down to `spacing`, without constructing the finest full
volume. Coil and component controls fetch only the selected data. Julia must remain
running; the default plot is self-contained.

# Keywords — both modes
All keywords in the signature apply to both modes. `spacing` is the embedded grid
spacing with `adaptive=false`, and the finest grid spacing with `adaptive=true`.
`fov`, `height`, `width`, and `darkmode` retain the same meaning in either mode.
"""
function plot_coil_sens(sys; fov=nothing, spacing=0.01, height=700, width=nothing,
    darkmode=false, adaptive=false)
    adaptive && return coil_spatial_plot(sys; fov, spacing, height, width, darkmode)
    receiver = sys.receiver
    fov = isnothing(fov) ? get_coil_sens_fov(receiver) : fov
    bounds = map(coil_plot_bounds(receiver), fov) do (lower, upper), extent
        center = (lower + upper) / 2
        (center - extent / 2, center + extent / 2)
    end
    grid = map(bounds, fov) do (lower, upper), extent
        intervals = floor(Int, extent / spacing)
        start = (lower + upper - intervals * spacing) / 2
        range(start; step=spacing, length=intervals + 1)
    end
    positions = [position for position in Iterators.product(grid...)
        if inside_coil_plot(receiver, position)]
    x, y, z = (getindex.(positions, dimension) for dimension in 1:3)
    sensitivities = get_sens(receiver, x, y, z)
    magnitude = [abs.(coil) for coil in eachcol(sensitivities)]
    phase = [angle.(coil) for coil in eachcol(sensitivities)]
    coordinates = (100 .* x, 100 .* y, 100 .* z)
    ncoils = get_n_coils(receiver)
    dimensions = findall(bound -> first(bound) != last(bound), bounds)
    view_2d = length(dimensions) < 3
    displayed = [dimensions; setdiff(1:3, dimensions)][1:2]
    axis_names = ("x", "y", "z")
    background, text, plot_background, grid_color, _ = theme_chooser(darkmode)
    magnitude_scale = attr(; cmin=0, cmax=maximum(abs, sensitivities), colorscale="Greys",
        colorbar=attr(; title=attr(; text="|S|"), thickness=16))
    phase_scale = attr(; cmin=-π, cmax=π, colorscale=PHASE_COLORSCALE,
        colorbar=attr(; title=attr(; text="rad"), thickness=16,
            tickvals=[-π, -π/2, 0, π/2, π], ticktext=["−π", "−π/2", "0", "π/2", "π"]))
    hover = view_2d ?
        "$(axis_names[displayed[1]]): %{x:.1f} cm<br>$(axis_names[displayed[2]]): %{y:.1f} cm<br>" :
        "x: %{x:.1f} cm<br>y: %{y:.1f} cm<br>z: %{z:.1f} cm<br>"
    traces = [
        (view_2d ? scattergl : scatter3d)(;
            x=coordinates[view_2d ? displayed[1] : 1],
            y=coordinates[view_2d ? displayed[2] : 2],
            (view_2d ? (;) : (; z=coordinates[3]))...,
            mode="markers", visible=coil == 1, showlegend=false, name="Coil $coil",
            marker=attr(; size=view_2d ? 4 : 2, color=magnitude[coil], coloraxis="coloraxis"),
            hovertemplate=hover * "|S|: %{marker.color:.3g}<extra>%{fullData.name}</extra>",
        ) for coil in 1:ncoils
    ]
    extent = max(spacing * 100, maximum(bound -> maximum(abs, bound), bounds) * 100)
    axes = [attr(; title=name, range=[-extent, extent], ticksuffix=" cm",
        backgroundcolor=plot_background, gridcolor=grid_color, zerolinecolor=grid_color)
        for name in axis_names]
    layout = Layout(;
        paper_bgcolor=background, plot_bgcolor=plot_background, font_color=text,
        margin=attr(; l=0, r=80, t=ncoils > 1 ? 100 : 30, b=90),
        coloraxis=magnitude_scale,
        uirevision="coil-sensitivities",
        sliders=ncoils == 1 ? [] : [attr(;
            active=0, x=0.12, len=0.78, y=1.04, yanchor="bottom",
            currentvalue=attr(; prefix="Coil: "),
            steps=[attr(; label=string(coil), method="restyle",
                args=[attr(; visible=collect(1:ncoils) .== coil)])
                for coil in 1:ncoils],
        )],
        updatemenus=[attr(;
            name="coil-component", type="buttons", direction="right",
            x=0.5, xanchor="center", y=0, yanchor="top", pad=attr(; t=30),
            bgcolor="#2a7fb8", bordercolor="#2a7fb8", font=attr(; color="#111111"),
            showactive=true, active=0,
            buttons=[attr(; label, method="update", args=[
                Dict("marker.color" => values, "hovertemplate" => hover * value_label *
                    ": %{marker.color:.3g}<extra>%{fullData.name}</extra>"),
                attr(; coloraxis=scale),
            ]) for (label, value_label, values, scale) in (
                ("Magnitude", "|S|", magnitude, magnitude_scale),
                ("Phase", "Phase (rad)", phase, phase_scale),
            )],
        )],
        modebar=attr(; orientation="h", bgcolor=background, color=text, activecolor=plot_background),
    )
    if view_2d
        layout[:xaxis] = axes[displayed[1]]
        layout[:yaxis] = axes[displayed[2]]
        layout[:xaxis][:scaleanchor] = "y"
    else
        layout[:scene] = attr(; xaxis=axes[1], yaxis=axes[2], zaxis=axes[3],
            aspectmode="manual", aspectratio=attr(; x=1, y=1, z=1))
    end
    height !== nothing && (layout.height = height)
    width !== nothing && (layout.width = width)
    config = PlotConfig(; displaylogo=false, toImageButtonOptions=attr(; format="svg").fields,
        modeBarButtonsToRemove=["zoom", "pan", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"])
    return PlotlyBase.Plot(traces, layout; config)
end

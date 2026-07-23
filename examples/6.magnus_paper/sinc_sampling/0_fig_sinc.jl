using CairoMakie

const OUTPUT_DIR = joinpath(@__DIR__, "output")
const FIGURE_FILE = joinpath(OUTPUT_DIR, "sinc_sampling.svg")
const ZOOM_CONSTANT_FILE = joinpath(OUTPUT_DIR, "sinc_sampling_zoom_constant.svg")
const ZOOM_LINEAR_FILE = joinpath(OUTPUT_DIR, "sinc_sampling_zoom_linear.svg")
const ZOOM_GL_FILE = joinpath(OUTPUT_DIR, "sinc_sampling_zoom_gl.svg")

const Δt = 0.6
const t0 = Δt
const t1 = 2Δt
const c_minus = 1 / 2 - sqrt(3) / 6
const c_plus = 1 / 2 + sqrt(3) / 6
const t_gl = [t0 + c_minus * Δt, t0 + c_plus * Δt]

sinc_rf(t) = sinc(1.1t)
linear_interp(t, t0, y0, t1, y1) = y0 + (y1 - y0) * (t - t0) / (t1 - t0)

function interval_arrow!(ax, y; shaftwidth=9, tipwidth=24, tiplength=24)
    arrows2d!(
        ax,
        [t0],
        [y],
        [Δt],
        [0.0];
        shaftwidth,
        tipwidth,
        tiplength,
        tail=Point2f[(0, 0), (1, -0.5), (1, 0.5)],
        taillength=tiplength,
        tailwidth=tipwidth,
        color=:black,
    )
    return nothing
end

function samples()
    return (;
        y0=sinc_rf(t0),
        y1=sinc_rf(t1),
        y_gl=sinc_rf.(t_gl),
        colors=Makie.wong_colors(),
    )
end

function render_figure()
    ts = range(-1.8, 1.58; length=1200)
    interval_t = range(t0, t1; length=200)
    (; y0, y1, y_gl, colors) = samples()

    fig = Figure(size=(920, 430), fontsize=42, figure_padding=6)
    ax = Axis(fig[1, 1])

    lines!(ax, ts, sinc_rf.(ts); linewidth=10, color=:black)
    hlines!(ax, [0.0]; linewidth=8, color=(:gray70, 0.35))
    foreach(t -> lines!(ax, [t, t], [-0.32, 1.06]; linewidth=4, color=(:gray45, 0.5)), (t0, t1))

    lines!(ax, [t0, t1], [y0, y0]; linewidth=9, color=colors[2])
    lines!(ax, interval_t, linear_interp.(interval_t, t0, y0, t1, y1); linewidth=9, color=colors[3])
    interval_arrow!(ax, -0.27; shaftwidth=9, tipwidth=28, tiplength=28)

    scatter!(ax, [t0], [y0]; markersize=48, color=colors[2])
    scatter!(ax, [t1], [y1]; markersize=48, color=colors[3])
    scatter!(ax, t_gl, y_gl; markersize=48, color=colors[4])
    xlims!(ax, -1.95, 1.72)
    ylims!(ax, -0.38, 1.12)
    hidedecorations!(ax)
    hidespines!(ax)
    resize_to_layout!(fig)
    return fig
end

function draw_zoom_base()
    ts = range(t0 - 0.35Δt, t1 + 0.20Δt; length=600)
    fig = Figure(size=(420, 430), fontsize=42, figure_padding=20)
    ax = Axis(fig[1, 1])

    lines!(ax, ts, sinc_rf.(ts); linewidth=8, color=:black)
    hlines!(ax, [0.0]; linewidth=6, color=(:gray70, 0.35))
    foreach(t -> lines!(ax, [t, t], [-0.30, 0.42]; linewidth=3, color=(:gray45, 0.5)), (t0, t1))

    interval_arrow!(ax, -0.25)
    xlims!(ax, t0 - 0.22Δt, t1 + 0.22Δt)
    ylims!(ax, -0.40, 0.50)
    hidedecorations!(ax)
    hidespines!(ax)
    return fig, ax
end

function render_zoom_constant()
    (; y0, colors) = samples()
    fig, ax = draw_zoom_base()
    lines!(ax, [t0, t1], [y0, y0]; linewidth=9, color=colors[2])
    scatter!(ax, [t0], [y0]; markersize=34, color=colors[2])
    resize_to_layout!(fig)
    return fig
end

function render_zoom_linear()
    interval_t = range(t0, t1; length=200)
    (; y0, y1, colors) = samples()
    fig, ax = draw_zoom_base()
    lines!(ax, interval_t, linear_interp.(interval_t, t0, y0, t1, y1); linewidth=9, color=colors[3])
    scatter!(ax, [t0, t1], [y0, y1]; markersize=34, color=colors[3])
    resize_to_layout!(fig)
    return fig
end

function render_zoom_gl()
    (; y_gl, colors) = samples()
    fig, ax = draw_zoom_base()
    scatter!(ax, t_gl, y_gl; markersize=34, color=colors[4])
    resize_to_layout!(fig)
    return fig
end

function (@main)(args)
    mkpath(OUTPUT_DIR)
    save(FIGURE_FILE, render_figure());
    save(ZOOM_CONSTANT_FILE, render_zoom_constant());
    save(ZOOM_LINEAR_FILE, render_zoom_linear());
    save(ZOOM_GL_FILE, render_zoom_gl());
    println(FIGURE_FILE)
    println(ZOOM_CONSTANT_FILE)
    println(ZOOM_LINEAR_FILE)
    println(ZOOM_GL_FILE)
end

using CairoMakie
using JLD2
using KomaMRI
using LaTeXStrings
using Statistics

include(joinpath(@__DIR__, "0_acc_bench.jl"))
const AB = AccBenchmarkCommon

const RESULT_FILE = joinpath(AB.OUTPUT_DIR, "convergence_all_methods_results.jld2")
const FIGURE_FILE = joinpath(AB.OUTPUT_DIR, "fig_convergence_all_methods.svg")
const FIGURE_PNG_FILE = joinpath(AB.OUTPUT_DIR, "fig_convergence_all_methods.png")
const SLOPE_FIT_POINTS = 7

all_magnus_methods() = (
    ("Const1", BlochMagnusConst1()),
    ("Lin2", BlochMagnusLin2()),
    ("Mid2", BlochMagnusMid2()),
    ("LinComm2", BlochMagnusLinComm2()),
    ("Quad2", BlochMagnusQuad2()),
    ("Quad4", BlochMagnusQuad4()),
    ("GL2", BlochMagnusGL2()),
    ("GL4", BlochMagnusGL4()),
    ("BGL4", BlochMagnusBGL4()),
    ("BGL6", BlochMagnusBGL6()),
)

function convergence_series(methods, obj, seq, sys, reference; precision)
    err = zeros(length(methods), length(AB.DT_RF))
    npoints = length(methods) * length(AB.DT_RF)
    for (method_i, (label, method)) in enumerate(methods)
        for (dt_i, Δt_rf) in enumerate(AB.DT_RF)
            point_i = (method_i - 1) * length(AB.DT_RF) + dt_i
            @info "convergence" point="$point_i/$npoints" precision method=label Δt_rf_us=Δt_rf * 1e6
            params = AB.sim_params(method, Δt_rf; gpu=false, precision)
            state = simulate(obj, seq, sys; sim_params=copy(params), verbose=false)
            err[method_i, dt_i] = AB.nrmse(state, reference)
        end
    end
    return err
end

function convergence_result()
    mkpath(AB.OUTPUT_DIR)
    seq = AB.hs_sequence()
    obj = AB.phantom_grid()
    sys = Scanner()
    methods = all_magnus_methods()
    @info "reference" precision="f64" method="BGL6" Δt_rf_us=AB.REFERENCE_DT_RF * 1e6
    reference = AB.reference_state(obj, seq, sys)
    return Dict(
        "duration" => AB.HS_DURATION,
        "native_samples" => AB.HS_NATIVE_SAMPLES,
        "nspins" => AB.NSPINS,
        "reference_dt_rf" => AB.REFERENCE_DT_RF,
        "dt_rf_us" => 1e6 .* AB.DT_RF,
        "method_labels" => collect(first.(methods)),
        "err64" => convergence_series(methods, obj, seq, sys, reference; precision="f64"),
        "err32" => convergence_series(methods, obj, seq, sys, reference; precision="f32"),
    )
end

function convergence_slope(dt_rf_us, err)
    inds = findall(>(0), err)[1:SLOPE_FIT_POINTS]
    x = log10.(dt_rf_us[inds])
    y = log10.(err[inds])
    return sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)).^2)
end

function render_figure(result)
    labels = result["method_labels"]
    dt_rf_us = result["dt_rf_us"]
    err64 = result["err64"]
    err32 = result["err32"]
    colors = [
        RGBf(0.88, 0.47, 0.10),
        RGBf(0.95, 0.65, 0.10),
        RGBf(0.55, 0.75, 0.35),
        RGBf(0.20, 0.60, 0.30),
        RGBf(0.25, 0.50, 0.85),
        RGBf(0.10, 0.25, 0.70),
        RGBf(0.85, 0.35, 0.70),
        RGBf(0.60, 0.20, 0.65),
        RGBf(0.55, 0.35, 0.20),
        RGBf(0.05, 0.05, 0.05),
    ]

    theme = Makie.theme_latexfonts()
    theme.fontsize = 30
    theme.linewidth = 3
    theme.markersize = 12

    with_theme(theme) do
        fig = Figure(size=(1050, 1220), figure_padding=(10, 10, 10, 10))
        ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, xreversed=true,
            xlabel=L"\Delta t\;(\mu\mathrm{s})", ylabel="NRMSE",
            xticks=([300, 100, 50, 20, 10, 5, 2, 1], string.([300, 100, 50, 20, 10, 5, 2, 1])),
            title="HS convergence: all methods", aspect=1)

        method_handles = LineElement[]
        for (i, (label, color)) in enumerate(zip(labels, colors))
            inds64 = findall(>(0), view(err64, i, :))
            inds32 = findall(>(0), view(err32, i, :))
            scatterlines!(ax, dt_rf_us[inds64], err64[i, inds64]; color)
            scatterlines!(ax, dt_rf_us[inds32], err32[i, inds32]; color, linestyle=:dash)
            push!(method_handles, LineElement(; color, linewidth=4))
        end

        ax.xminorticksvisible = true
        ax.yminorticksvisible = true
        ax.xminorticks = IntervalsBetween(10)
        ax.yminorticks = IntervalsBetween(10)
        ylims!(ax, 1e-13, 4e-1)
        xlims!(ax, maximum(dt_rf_us) * 1.15, minimum(dt_rf_us) / 1.15)

        Legend(fig[2, 1], method_handles, labels;
            orientation=:horizontal, tellwidth=false, nbanks=2)
        Legend(fig[3, 1],
            [LineElement(color=:black, linestyle=:solid, linewidth=4),
             LineElement(color=:black, linestyle=:dash, linewidth=4)],
            ["Float64", "Float32"];
            orientation=:horizontal, tellwidth=false, patchsize=(90, 30))
        slopes = [convergence_slope(dt_rf_us, view(err64, i, :)) for i in eachindex(labels)]
        slope_grid = fig[4, 1] = GridLayout()
        Label(slope_grid[1, 1:length(labels)], "Float64 fitted slope, Δt = 300-20 μs";
            fontsize=24, font=:bold, tellheight=true)
        for (i, (label, color, slope)) in enumerate(zip(labels, colors, slopes))
            Label(slope_grid[2, i], label; color, fontsize=21, tellheight=true)
            Label(slope_grid[3, i], string(round(slope; digits=2)); fontsize=21, tellheight=true)
        end
        rowsize!(fig.layout, 2, Fixed(95))
        rowsize!(fig.layout, 3, Fixed(70))
        rowsize!(fig.layout, 4, Fixed(105))
        return fig
    end
end

function main()
    result = isfile(RESULT_FILE) ? load(RESULT_FILE, "result") : convergence_result()
    isfile(RESULT_FILE) || jldsave(RESULT_FILE; result)
    fig = render_figure(result)
    save(FIGURE_FILE, fig)
    save(FIGURE_PNG_FILE, fig, px_per_unit=2)
    println(FIGURE_FILE)
    println(FIGURE_PNG_FILE)
    return FIGURE_FILE
end

main()

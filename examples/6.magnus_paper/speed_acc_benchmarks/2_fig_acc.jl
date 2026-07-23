# Figure: HS pulse speed/accuracy comparison for Bloch-Magnus methods.
#
# Run from the KomaMRI.jl root after generating one or both result files:
#   julia --project=examples/6.magnus_paper examples/6.magnus_paper/speed_acc_benchmarks/2_fig_acc.jl

using CairoMakie
using JLD2: load
using LaTeXStrings
using Printf: @sprintf
using Statistics

const OUTPUT_DIR = joinpath(@__DIR__, "output")
const CONVERGENCE_RESULTS_FILE = joinpath(OUTPUT_DIR, "convergence_acc_results.jld2")
const METAL_RESULTS_FILE = joinpath(OUTPUT_DIR, "metal_acc_results.jld2")
const CUDA_RESULTS_FILE = joinpath(OUTPUT_DIR, "cuda_acc_results.jld2")
const FIGURE_FILE = joinpath(OUTPUT_DIR, "fig_acc.svg")
const PLOT_DT_TICKS_US = [300, 100, 50, 20, 10, 5, 2, 1]
const CONVERGENCE_YLIMS = (1e-13, 4e-1)
const TIME_TICKS = (
    [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
    ["0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2", "5"],
)
const METHOD_SPEEDUP_TICKS = (
    10.0 .^ (0:3:9),
    ["1", rich("10", superscript("3")), rich("10", superscript("6")), rich("10", superscript("9"))],
)

errors(series, i) = view(series["err"], i, :)
time_samples_s(series, i) = collect(view(series["simulate_s_samples"], i, :))

order_label(order) = order == 1 ? L"\mathcal{O}(\Delta t)" : latexstring("\\mathcal{O}(\\Delta t^{$order})")

function order_label_geometry(dt_rf_us, err, label_dt_us, i1, i2; offset=0.045)
    xspan = log10(maximum(dt_rf_us) / minimum(dt_rf_us))
    yspan = log10(CONVERGENCE_YLIMS[2] / CONVERGENCE_YLIMS[1])
    u1 = log10(maximum(dt_rf_us) / dt_rf_us[i1]) / xspan
    u2 = log10(maximum(dt_rf_us) / dt_rf_us[i2]) / xspan
    v1 = log10(err[i1] / CONVERGENCE_YLIMS[1]) / yspan
    v2 = log10(err[i2] / CONVERGENCE_YLIMS[1]) / yspan
    m = (v2 - v1) / (u2 - u1)
    u = log10(maximum(dt_rf_us) / label_dt_us) / xspan
    v = v1 + m * (u - u1)
    du, dv = (-m, 1.0) ./ hypot(m, 1.0) .* offset
    x = maximum(dt_rf_us) / 10.0^((u + du) * xspan)
    y = CONVERGENCE_YLIMS[1] * 10.0^((v + dv) * yspan)
    return x, y, atan(m)
end

function pre_roundoff_indices(err)
    i = something(findfirst(i -> err[i + 1] >= err[i], 1:(length(err) - 1)), length(err))
    return (i > 2 && i < length(err)) ? (1:i-1) : (1:i)
end

function spread(xs)
    q25, q50, q75 = quantile(xs, (0.25, 0.50, 0.75))
    return q50, q50 - q25, q75 - q50
end
spreads(xs) = getindex.(spread.(xs), 1), getindex.(spread.(xs), 2), getindex.(spread.(xs), 3)
time_spreads_s(series, i) = spreads(time_samples_s(series, i))

function gpu_speedup_at_same_step(err, cpu_samples, gpu_samples)
    inds = pre_roundoff_indices(err)
    err = err[inds]
    cpu_samples = cpu_samples[inds]
    gpu_samples = gpu_samples[inds]
    speedup = zeros(length(err))
    low = similar(speedup)
    high = similar(speedup)
    for i in eachindex(err)
        speedups = [cpu / gpu for cpu in cpu_samples[i] for gpu in gpu_samples[i]]
        speedup[i], low[i], high[i] = spread(speedups)
    end
    return err, speedup, low, high
end

function baseline_model(err, samples)
    inds = pre_roundoff_indices(err)
    err = collect(err[inds])
    samples = samples[inds]
    order = sortperm(err)
    return err[order], samples[order]
end

function extrapolated_time_samples(err, baseline_err, baseline_samples)
    i = clamp(searchsortedlast(baseline_err, err), 1, length(baseline_err) - 1)
    w = (log10(err) - log10(baseline_err[i])) / (log10(baseline_err[i + 1]) - log10(baseline_err[i]))
    return [
        10.0^((1 - w) * log10(t1) + w * log10(t2))
        for t1 in baseline_samples[i]
        for t2 in baseline_samples[i + 1]
    ]
end

function method_speedup_over_baseline(err, samples, baseline_err, baseline_samples)
    inds = pre_roundoff_indices(err)
    err = collect(err[inds])
    samples = samples[inds]
    baseline_err, baseline_samples = baseline_model(baseline_err, baseline_samples)
    speedup = zeros(length(err))
    low = similar(speedup)
    high = similar(speedup)
    for i in eachindex(err)
        values = [
            baseline_time / t
            for baseline_time in extrapolated_time_samples(err[i], baseline_err, baseline_samples)
            for t in samples[i]
        ]
        speedup[i], low[i], high[i] = spread(values)
    end
    return err, speedup, low, high
end

nrmse_ticks(kmax=7; step=1) = (10.0 .^ (-1:-step:-kmax), [rich("10", superscript("-$k")) for k in 1:step:kmax])
const NRMSE_TICKS = nrmse_ticks(13; step=3)

function tick_label(x)
    x < 10 && return @sprintf("%.2g", x)
    x < 1000 && return @sprintf("%.0f", x)
    return replace(@sprintf("%.1e", x), "e+0" => "e", "e+" => "e", "e-0" => "e-")
end

function speedup_ticks(values; min_log_gap=0.0)
    ticks = sort(unique(filter(x -> isfinite(x) && x > 0, values)))
    if min_log_gap > 0
        spaced_ticks = Float64[]
        for tick in ticks
            if isempty(spaced_ticks) || log10(tick) - log10(last(spaced_ticks)) >= min_log_gap
                push!(spaced_ticks, tick)
            else
                spaced_ticks[end] = tick
            end
        end
        ticks = spaced_ticks
    end
    return ticks, tick_label.(ticks)
end

function gpu_speedup_ticks(series, method_indices; extra_values=Float64[])
    values = [1.0; extra_values]
    all_values = [y for s in series for y in s.speedup]
    isempty(all_values) || push!(values, minimum(all_values))
    for i in method_indices
        vals = [y for s in series if s.method == i for y in s.speedup]
        isempty(vals) || push!(values, maximum(vals))
    end
    return speedup_ticks(values; min_log_gap=0.05)
end

function speedup_hi(series)
    return 1.2 * maximum((
        y for s in series for y in s.speedup .+ s.high
        if isfinite(y) && y > 0
    ); init=1.0)
end

function gpu_speedup_series(result, method_indices, colors)
    result === nothing && return NamedTuple[]
    series = NamedTuple[]
    cpu32 = result["cpu32"]
    cpu64 = result["cpu64"]
    gpu32 = result["gpu32"]
    gpu64 = result["gpu64"]
    for (i, color) in zip(method_indices, colors)
        err, speedup, low, high = gpu_speedup_at_same_step(errors(cpu32, i), time_samples_s(cpu32, i), time_samples_s(gpu32, i))
        push!(series, (; method=i, err, speedup, low, high, color, linestyle=:dash))
        if gpu64 !== nothing
            err, speedup, low, high = gpu_speedup_at_same_step(errors(cpu64, i), time_samples_s(cpu64, i), time_samples_s(gpu64, i))
            push!(series, (; method=i, err, speedup, low, high, color, linestyle=:solid))
        end
    end
    return series
end

function time_series(result, fields, method_indices, colors)
    result === nothing && return NamedTuple[]
    series = NamedTuple[]
    for (field, linestyle) in fields
        data = result[field]
        data === nothing && continue
        for (i, color) in zip(method_indices, colors)
            method_err = errors(data, i)
            inds = pre_roundoff_indices(method_err)
            time, low, high = getindex.(time_spreads_s(data, i), Ref(inds))
            push!(series, (; method=i, err=method_err[inds], time, low, high, color, linestyle))
        end
    end
    return series
end

function method_speedup_series(cpu64, cpu32, method_indices, colors)
    series = NamedTuple[]
    baseline_i = first(method_indices)
    for (cpu, linestyle) in ((cpu64, :solid), (cpu32, :dash))
        baseline_err = errors(cpu, baseline_i)
        baseline_samples = time_samples_s(cpu, baseline_i)
        for (i, color) in zip(method_indices[2:end], colors[2:end])
            err, speedup, low, high = method_speedup_over_baseline(
                errors(cpu, i),
                time_samples_s(cpu, i),
                baseline_err,
                baseline_samples,
            )
            push!(series, (; method=i, err, speedup, low, high, color, linestyle))
        end
    end
    return series
end

function convergence_accuracy_limits(conv_err64, conv_err32, method_indices, colors)
    limits = NamedTuple[]
    for (i, color) in zip(method_indices, colors)
        for (err, linestyle) in ((view(conv_err64, i, :), :solid), (view(conv_err32, i, :), :dash))
            positive_err = [x for x in err if isfinite(x) && x > 0]
            isempty(positive_err) || push!(limits, (; method=i, err=minimum(positive_err), color, linestyle))
        end
    end
    return limits
end

function baseline_extrapolation_series(cpu64, cpu32, method_indices, colors)
    series = NamedTuple[]
    baseline_i = first(method_indices)
    for cpu in (cpu64, cpu32)
        baseline_err, baseline_samples = baseline_model(errors(cpu, baseline_i), time_samples_s(cpu, baseline_i))
        target_err = minimum((
            minimum(errors(cpu, i)[pre_roundoff_indices(errors(cpu, i))])
            for i in method_indices[2:end]
        ); init=minimum(baseline_err))
        target_err < minimum(baseline_err) || continue
        err = 10.0 .^ range(log10(minimum(baseline_err)), log10(target_err); length=80)
        time = first.(spread.(extrapolated_time_samples.(err, Ref(baseline_err), Ref(baseline_samples))))
        push!(series, (; err, time, color=first(colors)))
    end
    return series
end

function time_limits(series)
    values = [y for s in series for y in s.time if isfinite(y) && y > 0]
    isempty(values) && return (1e-3, 1.0)
    return minimum(values) / 1.15, maximum(values) * 1.15
end

function bandlines!(ax, x, y, low, high; color, linestyle=:solid, marker=:circle)
    inds = findall(i -> all(isfinite, (x[i], y[i], low[i], high[i])) && y[i] > 0 && low[i] >= 0 && high[i] >= 0, eachindex(x))
    isempty(inds) && return nothing
    band_inds = [i for i in inds if low[i] > 0 || high[i] > 0]
    if !isempty(band_inds)
        ylo = max.(y[band_inds] .- low[band_inds], floatmin(Float64))
        yhi = y[band_inds] .+ high[band_inds]
        band!(ax, x[band_inds], ylo, yhi; color=(color, 0.18))
    end
    scatterlines!(ax, x[inds], y[inds]; color, linestyle, marker)
    return nothing
end

function plot_time_series!(ax, series)
    for s in series
        bandlines!(ax, s.err, s.time, s.low, s.high; color=s.color, linestyle=s.linestyle)
    end
    return nothing
end

function render_figure(convergence, metal, cuda)
    data = metal !== nothing ? metal : cuda
    data === nothing && error("Generate metal_acc_results.jld2 or cuda_acc_results.jld2 first.")
    convergence === nothing && error("Generate convergence_acc_results.jld2 first.")
    for result in (metal, cuda), key in ("duration", "native_samples", "nspins", "reference_dt_rf", "dt_rf_us", "method_labels")
        result === nothing || convergence[key] == result[key] || error("result files use different $key")
    end

    cpu64 = data["cpu64"]
    cpu32 = data["cpu32"]
    labels = data["method_labels"]
    method_indices = eachindex(labels)
    colors = Makie.wong_colors()[[2, 3, 4, 6]]
    conv_dt_rf_us = convergence["dt_rf_us"]
    conv_err64 = convergence["err64"]
    conv_err32 = convergence["err32"]
    metal_gpu_speedups = gpu_speedup_series(metal, method_indices, colors)
    cuda_gpu_speedups = gpu_speedup_series(cuda, method_indices, colors)
    cpu_times = time_series(data, (("cpu64", :solid), ("cpu32", :dash)), method_indices, colors)
    baseline_extrapolations = baseline_extrapolation_series(cpu64, cpu32, method_indices, colors)
    method_speedups = method_speedup_series(cpu64, cpu32, method_indices, colors)
    convergence_accuracy = convergence_accuracy_limits(conv_err64, conv_err32, method_indices, colors)
    simulation_time_lims = time_limits(cpu_times)

    metal_speedup_hi = speedup_hi(metal_gpu_speedups)
    cuda_speedup_hi = speedup_hi(cuda_gpu_speedups)
    gpu_speedup_hi = max(metal_speedup_hi, cuda_speedup_hi)
    cuda_m46_f32_speedup = maximum((
        maximum(s.speedup) for s in cuda_gpu_speedups
        if s.linestyle == :dash && (s.method == 3 || s.method == 4)
    ); init=-Inf)

    theme = Makie.theme_latexfonts()
    theme.fontsize = 32
    theme.labelsize = 32
    theme.linewidth = 4
    theme.markersize = 18

    with_theme(theme) do
        fig = Figure(size=(1800, 1900), figure_padding=(0, 35, 0, 0))
        ax_a = Axis(fig[1, 1]; xscale=log10, yscale=log10, xreversed=true,
            xlabel=rich("Δ", rich("t"; font=:italic), " (μs)"), ylabel="NRMSE",
            xticks=(PLOT_DT_TICKS_US, string.(PLOT_DT_TICKS_US)), yticks=NRMSE_TICKS,
            aspect=1, title="(A) Method convergence")
        ax_b = Axis(fig[1, 2]; xscale=log10, yscale=log10, xreversed=true,
            xlabel="NRMSE", ylabel="Simulation time (s)", xticks=NRMSE_TICKS,
            yticks=TIME_TICKS, aspect=1, title="(B) CPU accuracy-cost tradeoff")

        ax_c = Axis(fig[2, 1]; xscale=log10, yscale=log10, xreversed=true,
            xlabel="NRMSE", ylabel="Speedup over Magnus1",
            xticks=NRMSE_TICKS, yticks=METHOD_SPEEDUP_TICKS,
            yticklabelsize=22, aspect=1, title="(C) Method speedup")

        gpu_grid = fig[2, 2] = GridLayout()
        Label(fig[2, 2, Top()], "(D) GPU speedup over CPU"; fontsize=32, font=:bold, tellwidth=false)
        Label(gpu_grid[1, 1], metal === nothing ? "Metal (Apple Silicon)" : metal["gpu_title"]; fontsize=30, font=:bold, tellwidth=false)
        Label(gpu_grid[1, 2], "NVIDIA A100"; fontsize=30, font=:bold, tellwidth=false)
        ax_d_metal = Axis(gpu_grid[2, 1]; xscale=log10, yscale=log10, xreversed=true,
            xlabel="NRMSE", ylabel="Speedup", xticks=NRMSE_TICKS, xticklabelsize=26,
            yticks=gpu_speedup_ticks(metal_gpu_speedups, method_indices),
            yticklabelsize=20)
        ax_d_cuda = Axis(gpu_grid[2, 2]; xscale=log10, yscale=log10, xreversed=true,
            xlabel="NRMSE", xticks=NRMSE_TICKS, xticklabelsize=26,
            yticks=gpu_speedup_ticks(cuda_gpu_speedups, method_indices; extra_values=[cuda_m46_f32_speedup]),
            yticklabelsize=20)
        colgap!(gpu_grid, 18)
        rowsize!(gpu_grid, 1, Fixed(38))

        for (i, label, color) in zip(method_indices, labels, colors)
            method_err64 = view(conv_err64, i, :)
            method_err32 = view(conv_err32, i, :)
            inds64 = findall(>(0), method_err64)
            inds32 = findall(>(0), method_err32)
            scatterlines!(ax_a, conv_dt_rf_us[inds64], method_err64[inds64]; color, label)
            scatterlines!(ax_a, conv_dt_rf_us[inds32], method_err32[inds32]; color, linestyle=:dash)
        end
        label_dt_us = 35.0
        i50 = findfirst(==(50.0), conv_dt_rf_us)
        i20 = findfirst(==(20.0), conv_dt_rf_us)
        for (i, order) in ((1, 1), (2, 2), (3, 4), (4, 6))
            err64 = view(conv_err64, i, :)
            label_x, label_y, rotation = order_label_geometry(conv_dt_rf_us, err64, label_dt_us, i50, i20)
            text!(ax_a, label_x, label_y; text=order_label(order), align=(:center, :center),
                rotation, fontsize=32, color=:black)
        end

        plot_time_series!(ax_b, cpu_times)
        for s in baseline_extrapolations
            lines!(ax_b, s.err, s.time; color=s.color, linestyle=:dot, linewidth=4)
        end
        for s in method_speedups
            bandlines!(ax_c, s.err, s.speedup, s.low, s.high; color=s.color, linestyle=s.linestyle)
        end
        for s in convergence_accuracy
            is_baseline = s.method == first(method_indices)
            vlines!(ax_c, [s.err]; color=(s.color, is_baseline ? 0.8 : 0.45),
                linestyle=s.linestyle, linewidth=is_baseline ? 4 : 2)
        end

        for (ax, series) in ((ax_d_metal, metal_gpu_speedups), (ax_d_cuda, cuda_gpu_speedups))
            for s in series
                bandlines!(ax, s.err, s.speedup, s.low, s.high; color=s.color, linestyle=s.linestyle)
            end
        end

        for ax in (ax_a, ax_b, ax_c, ax_d_metal, ax_d_cuda)
            ax.xminorticksvisible = true
            ax.yminorticksvisible = true
            ax.xminorticks = IntervalsBetween(10)
            ax.yminorticks = IntervalsBetween(10)
        end

        ylims!(ax_a, CONVERGENCE_YLIMS...)
        ylims!(ax_b, simulation_time_lims...)
        ylims!(ax_c, 1, 1e9)
        ylims!(ax_d_metal, 1, gpu_speedup_hi)
        ylims!(ax_d_cuda, 1, gpu_speedup_hi)
        xlims!(ax_b, 2e-1, 1e-13)
        xlims!(ax_c, 2e-1, 1e-13)
        xlims!(ax_d_metal, 2e-1, 1e-13)
        xlims!(ax_d_cuda, 2e-1, 1e-13)

        Legend(fig[3, 1:2], [LineElement(color=color, linewidth=4) for color in colors], collect(labels);
            orientation=:horizontal, patchsize=(70, 30), tellwidth=false, tellheight=true)
        precision_lines = [
            LineElement(color=:black, linestyle=:solid, linewidth=4),
            LineElement(color=:black, linestyle=:dash, linewidth=4),
            LineElement(color=colors[1], linestyle=:dot, linewidth=4),
        ]
        precision_labels = ["Float64", "Float32", "Magnus1 extrapolated"]
        Legend(fig[4, 1:2],
            precision_lines,
            precision_labels;
            orientation=:horizontal, patchsize=(70, 30), tellwidth=false, tellheight=true)
        rowsize!(fig.layout, 3, Fixed(65))
        rowsize!(fig.layout, 4, Fixed(65))
        return fig
    end
end

function (@main)(args)
    isempty(args) || error("expected no arguments")
    convergence = isfile(CONVERGENCE_RESULTS_FILE) ? load(CONVERGENCE_RESULTS_FILE, "result") : nothing
    metal = isfile(METAL_RESULTS_FILE) ? load(METAL_RESULTS_FILE, "result") : nothing
    cuda = isfile(CUDA_RESULTS_FILE) ? load(CUDA_RESULTS_FILE, "result") : nothing
    fig = render_figure(convergence, metal, cuda)
    mkpath(OUTPUT_DIR)
    save(FIGURE_FILE, fig);
    println(FIGURE_FILE)
end

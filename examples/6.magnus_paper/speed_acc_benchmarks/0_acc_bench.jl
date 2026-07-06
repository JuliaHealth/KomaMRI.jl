module AccBenchmarkCommon

using BenchmarkTools
using KomaMRI
using Statistics

export OUTPUT_DIR, benchmark_result, convergence_result, parse_benchmark_samples

const OUTPUT_DIR = joinpath(@__DIR__, "output")
const HS_DURATION = 18.3e-3
const HS_NATIVE_SAMPLES = 400_001
const NSPINS = 100_000
const DEFAULT_BENCHMARK_SAMPLES = 100
const BENCHMARK_SECONDS = 1e9
const REFERENCE_DT_RF = 1e-6
const DT_RF = [300e-6, 250e-6, 200e-6, 150e-6, 100e-6, 50e-6, 20e-6, 10e-6, 5e-6, 2e-6, 1e-6]

magnus_methods() = (
    ("Magnus1", BlochMagnus1()),
    ("Magnus2", BlochMagnus2()),
    ("Magnus4", BlochMagnus4()),
    ("Magnus6", BlochMagnus6()),
)

function hs_sequence(n=HS_NATIVE_SAMPLES)
    b1max = 13.5e-6
    β̂ = 4
    μ = 6
    β = 2β̂ / HS_DURATION
    t = range(-HS_DURATION / 2, HS_DURATION / 2, n)
    B1 = b1max .* sech.(β .* t)
    Δf = -μ * β .* tanh.(β .* t) ./ (2π)

    seq = Sequence()
    @addblock seq += RF(B1, HS_DURATION, Δf, 0)
    return seq
end

phantom_grid(nspins=NSPINS) = Phantom(
    x=collect(range(-6e-2, 6e-2; length=nspins)),
    y=zeros(nspins),
    z=zeros(nspins),
    ρ=ones(nspins),
    T1=fill(Inf, nspins),
    T2=fill(Inf, nspins),
    Δw=2π .* collect(range(-2e3, 2e3; length=nspins)),
)

sim_params(method, Δt_rf; gpu, precision) = Dict{String,Any}(
    "gpu" => gpu,
    "Nthreads" => gpu ? 1 : Threads.nthreads(),
    "return_type" => "state",
    "precision" => precision,
    "sim_method" => method,
    "Δt" => Inf,
    "Δt_rf" => Δt_rf,
)

function parse_benchmark_samples(args)
    length(args) <= 1 || error("expected at most one argument: benchmark_samples")
    isempty(args) && return DEFAULT_BENCHMARK_SAMPLES
    samples = parse(Int, only(args))
    samples > 0 || error("benchmark_samples must be positive")
    return samples
end

function nrmse(M, Mref)
    err2 = abs2.(M.xy .- Mref.xy) .+ abs2.(M.z .- Mref.z)
    ref2 = abs2.(Mref.xy) .+ abs2.(Mref.z)
    return sqrt(sum(err2) / sum(ref2))
end

function benchmark_simulate_s(obj, seq, sys, params, benchmark_samples)
    simulate(obj, seq, sys; sim_params=copy(params), verbose=false)
    trial = @benchmark simulate($obj, $seq, $sys; sim_params=copy($params), verbose=false) samples=benchmark_samples evals=1 seconds=BENCHMARK_SECONDS
    times = Float64.(trial.times)
    all(isfinite, times) || error("non-finite benchmark timing")
    length(times) == benchmark_samples || error("expected $benchmark_samples samples, got $(length(times))")
    times_s = times ./ 1e9
    return median(times_s), std(times_s), times_s
end

reference_state(obj, seq, sys; gpu=false) = simulate(
    obj,
    seq,
    sys;
    sim_params=sim_params(BlochMagnus6(), REFERENCE_DT_RF; gpu, precision="f64"),
    verbose=false,
)

function compute_series(obj, seq, sys, reference; gpu, precision, benchmark_samples)
    methods = magnus_methods()
    err = zeros(length(methods), length(DT_RF))
    simulate_s = zeros(length(methods), length(DT_RF))
    simulate_s_std = zeros(length(methods), length(DT_RF))
    simulate_s_samples = [Float64[] for _ in methods, _ in DT_RF]
    backend = gpu ? "GPU" : "CPU"
    npoints = length(methods) * length(DT_RF)
    for (method_i, (label, method)) in enumerate(methods)
        for (dt_i, Δt_rf) in enumerate(DT_RF)
            point_i = (method_i - 1) * length(DT_RF) + dt_i
            @info "benchmark" point="$point_i/$npoints" precision backend method=label Δt_rf_us=Δt_rf * 1e6
            t0 = time()
            params = sim_params(method, Δt_rf; gpu, precision)
            state = simulate(obj, seq, sys; sim_params=copy(params), verbose=false)
            err[method_i, dt_i] = nrmse(state, reference)
            simulate_s[method_i, dt_i], simulate_s_std[method_i, dt_i], simulate_s_samples[method_i, dt_i] =
                benchmark_simulate_s(obj, seq, sys, params, benchmark_samples)
            @info "benchmark done" point="$point_i/$npoints" elapsed_s=round(time() - t0; digits=1) median_s=simulate_s[method_i, dt_i] std_s=simulate_s_std[method_i, dt_i]
        end
    end
    return Dict(
        "err" => err,
        "simulate_s" => simulate_s,
        "simulate_s_std" => simulate_s_std,
        "simulate_s_samples" => simulate_s_samples,
    )
end

function convergence_series(obj, seq, sys, reference; precision)
    methods = magnus_methods()
    err = zeros(length(methods), length(DT_RF))
    npoints = length(methods) * length(DT_RF)
    for (method_i, (label, method)) in enumerate(methods)
        for (dt_i, Δt_rf) in enumerate(DT_RF)
            point_i = (method_i - 1) * length(DT_RF) + dt_i
            @info "convergence" point="$point_i/$npoints" precision method=label Δt_rf_us=Δt_rf * 1e6
            params = sim_params(method, Δt_rf; gpu=false, precision)
            state = simulate(obj, seq, sys; sim_params=copy(params), verbose=false)
            err[method_i, dt_i] = nrmse(state, reference)
        end
    end
    return err
end

function benchmark_result(; platform, gpu_title, include_gpu64, benchmark_samples=DEFAULT_BENCHMARK_SAMPLES, reference_gpu=false)
    seq = hs_sequence()
    obj = phantom_grid()
    sys = Scanner()
    @info "reference" backend=reference_gpu ? "GPU" : "CPU" precision="f64" method="Magnus6" Δt_rf_us=REFERENCE_DT_RF * 1e6
    t0 = time()
    reference = reference_state(obj, seq, sys; gpu=reference_gpu)
    @info "reference done" elapsed_s=round(time() - t0; digits=1)
    return Dict(
        "platform" => platform,
        "gpu_title" => gpu_title,
        "duration" => HS_DURATION,
        "native_samples" => HS_NATIVE_SAMPLES,
        "nspins" => NSPINS,
        "benchmark_samples" => benchmark_samples,
        "reference_dt_rf" => REFERENCE_DT_RF,
        "dt_rf_us" => 1e6 .* DT_RF,
        "method_labels" => collect(first.(magnus_methods())),
        "cpu64" => compute_series(obj, seq, sys, reference; gpu=false, precision="f64", benchmark_samples),
        "cpu32" => compute_series(obj, seq, sys, reference; gpu=false, precision="f32", benchmark_samples),
        "gpu32" => compute_series(obj, seq, sys, reference; gpu=true, precision="f32", benchmark_samples),
        "gpu64" => include_gpu64 ? compute_series(obj, seq, sys, reference; gpu=true, precision="f64", benchmark_samples) : nothing,
    )
end

function convergence_result(; nspins=NSPINS)
    seq = hs_sequence()
    obj = phantom_grid(nspins)
    sys = Scanner()
    reference = reference_state(obj, seq, sys)
    return Dict(
        "duration" => HS_DURATION,
        "native_samples" => HS_NATIVE_SAMPLES,
        "nspins" => nspins,
        "reference_dt_rf" => REFERENCE_DT_RF,
        "dt_rf_us" => 1e6 .* DT_RF,
        "method_labels" => collect(first.(magnus_methods())),
        "err64" => convergence_series(obj, seq, sys, reference; precision="f64"),
        "err32" => convergence_series(obj, seq, sys, reference; precision="f32"),
    )
end

end

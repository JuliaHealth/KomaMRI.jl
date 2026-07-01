# Benchmark HS pulse speed/accuracy on Metal (Apple Silicon).
#
# Run from the KomaMRI.jl root:
#   julia --project=examples/6.magnus_paper --threads=auto examples/6.magnus_paper/speed_acc_benchmarks/1_run_metal.jl

using Metal
using JLD2: jldsave

include("0_acc_bench.jl")
using .AccBenchmarkCommon

const RESULT_FILE = joinpath(OUTPUT_DIR, "metal_acc_results.jld2")

function (@main)(args)
    Metal.functional() || error("Metal is not functional")
    benchmark_samples = parse_benchmark_samples(args)
    result = benchmark_result(;
        platform=:metal,
        gpu_title=String(Metal.current_device().name),
        include_gpu64=false,
        benchmark_samples,
    )
    mkpath(OUTPUT_DIR)
    jldsave(RESULT_FILE; result);
    println(RESULT_FILE)
end

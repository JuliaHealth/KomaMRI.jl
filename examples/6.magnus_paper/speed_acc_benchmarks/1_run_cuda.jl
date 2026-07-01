# Benchmark HS pulse speed/accuracy on CUDA (NVIDIA GPU).
#
# Run from the KomaMRI.jl root:
#   julia --project=examples/6.magnus_paper --threads=auto examples/6.magnus_paper/speed_acc_benchmarks/1_run_cuda.jl

using CUDA
using JLD2: jldsave

include("0_acc_bench.jl")
using .AccBenchmarkCommon

const RESULT_FILE = joinpath(OUTPUT_DIR, "cuda_acc_results.jld2")

function (@main)(args)
    CUDA.functional() || error("CUDA is not functional")
    benchmark_samples = parse_benchmark_samples(args)
    result = benchmark_result(;
        platform=:cuda,
        gpu_title=CUDA.name(CUDA.device()),
        include_gpu64=true,
        benchmark_samples,
        reference_gpu=true,
    )
    mkpath(OUTPUT_DIR)
    jldsave(RESULT_FILE; result);
    println(RESULT_FILE)
end

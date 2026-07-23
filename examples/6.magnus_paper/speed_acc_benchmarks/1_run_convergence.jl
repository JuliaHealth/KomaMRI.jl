# Generate high-native-waveform data for panel A.

using JLD2: jldsave

include("0_acc_bench.jl")
using .AccBenchmarkCommon

const RESULT_FILE = joinpath(OUTPUT_DIR, "convergence_acc_results.jld2")

function (@main)(args)
    isempty(args) || error("expected no arguments")
    result = convergence_result()
    mkpath(OUTPUT_DIR)
    jldsave(RESULT_FILE; result);
    println(RESULT_FILE)
end

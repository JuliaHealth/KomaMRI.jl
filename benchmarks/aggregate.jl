using BenchmarkTools

const GPU_BACKENDS = ["AMDGPU", "CUDA", "Metal", "oneAPI"]
const NUM_CPU_THREADS = [1, 2, 4, 8]
const CPU_BENCHMARK_FILES = ["CPUbenchmarks$(n)threads.json" for n in NUM_CPU_THREADS]
const GPU_BENCHMARK_FILES = ["$(backend)benchmarks.json" for backend in GPU_BACKENDS]

function merge_results!(results::BenchmarkTools.BenchmarkGroup, new_results::BenchmarkTools.BenchmarkGroup)
    for key in keys(new_results)
        if haskey(results, key) &&
           results[key] isa BenchmarkTools.BenchmarkGroup &&
           new_results[key] isa BenchmarkTools.BenchmarkGroup
            merge_results!(results[key], new_results[key])
        else
            results[key] = new_results[key]
        end
    end
    return results
end

function load_benchmark_group(filepath)
    results = BenchmarkTools.load(filepath)[1]
    if results isa BenchmarkTools.BenchmarkGroup
        return results
    else
        @warn "Unexpected file format for file at path: $(filepath)"
        return nothing
    end
end

function aggregate_results(result_dir)
    results = BenchmarkTools.BenchmarkGroup()
    for filename in vcat(CPU_BENCHMARK_FILES, GPU_BENCHMARK_FILES)
        filepath = joinpath(result_dir, filename)
        if ispath(filepath)
            loaded_results = load_benchmark_group(filepath)
            if loaded_results !== nothing
                merge_results!(results, loaded_results)
            end
        else
            @warn "No file found at path: $(filepath)"
        end
    end
    @assert !isempty(results) "No benchmark results found"
    return results
end

const RESULT_DIR = joinpath(@__DIR__, "results")
const RESULTS = aggregate_results(RESULT_DIR)
BenchmarkTools.save(joinpath(RESULT_DIR, "combinedbenchmarks.json"), RESULTS)

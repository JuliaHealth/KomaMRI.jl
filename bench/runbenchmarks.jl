using KomaMRI
using Pkg
using Suppressor
using BenchmarkTools

const SUITE = BenchmarkGroup()

# To run benchmarks on a specific GPU backend, add AMDGPU / CUDA / Metal / oneAPI
# to bench/Project.toml and change the variable below to the backend name
const BENCHMARK_GROUP = get(ENV, "BENCHMARK_GROUP", "CPU")

# Number of CPU threads to benchmarks on 
const BENCHMARK_CPU_THREADS = (1,2,4,8)
if BENCHMARK_GROUP == "CPU" && Threads.nthreads() < maximum(BENCHMARK_CPU_THREADS)
    @error "More CPU threads were requested than are available. Change the JULIA_NUM_THREADS environment variable or pass --threads=$(maximum(BENCHMARK_CPU_THREADS)) as a julia argument"
end

if BENCHMARK_GROUP == "AMDGPU"
    using AMDGPU # ] add AMDGPU to bench/Project.toml 
    @info "Running AMDGPU benchmarks" maxlog=1
elseif BENCHMARK_GROUP == "CUDA"
    using CUDA # ] add CUDA to bench/Project.toml 
    @info "Running CUDA benchmarks" maxlog=1
elseif BENCHMARK_GROUP == "Metal"
    using Metal # ] add Metal to bench/Project.toml 
    @info "Running Metal benchmarks" maxlog=1
elseif BENCHMARK_GROUP == "oneAPI"
    using oneAPI # ] add oneAPI to bench/Project.toml 
    @info "Running oneAPI benchmarks" maxlog=1
else
    @info "Running CPU benchmarks with threads=$(BENCHMARK_CPU_THREADS)"
end

include("setup.jl")
setup_benchmarks(SUITE, BENCHMARK_GROUP, BENCHMARK_CPU_THREADS)
BenchmarkTools.tune!(SUITE; verbose=true)
results = BenchmarkTools.run(SUITE; verbose=true)
filepath = joinpath(dirname(@__FILE__), "results")
mkpath(filepath)
filename = string(BENCHMARK_GROUP, "benchmarks.json")
BenchmarkTools.save(joinpath(filepath, filename), median(results))
@info "Saved results to $(joinpath(filepath, filename))"
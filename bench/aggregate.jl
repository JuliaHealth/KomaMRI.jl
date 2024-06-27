using BenchmarkTools

const CPU_filepath = joinpath(dirname(@__FILE__), "results", "CPUbenchmarks.json")
const AMDGPU_filepath = joinpath(dirname(@__FILE__), "results", "AMDGPUbenchmarks.json")
const CUDA_filepath = joinpath(dirname(@__FILE__), "results", "CUDAbenchmarks.json")
const Metal_filepath = joinpath(dirname(@__FILE__), "results", "Metalbenchmarks.json")
const oneAPI_filepath = joinpath(dirname(@__FILE__), "results", "oneAPIbenchmarks.json")

@assert ispath(CPU_filepath)
@assert ispath(AMDGPU_filepath)
@assert ispath(CUDA_filepath)
@assert ispath(Metal_filepath)
@assert ispath(oneAPI_filepath)

const CPU_results = BenchmarkTools.load(CPU_filepath)[1]
const AMDGPU_results = BenchmarkTools.load(AMDGPU_filepath)[1]
const CUDA_results = BenchmarkTools.load(CUDA_filepath)[1]
const Metal_results = BenchmarkTools.load(Metal_filepath)[1]
const oneAPI_results = BenchmarkTools.load(oneAPI_filepath)[1]

@assert CPU_results isa BenchmarkTools.BenchmarkGroup
@assert AMDGPU_results isa BenchmarkTools.BenchmarkGroup
@assert CUDA_results isa BenchmarkTools.BenchmarkGroup
@assert Metal_results isa BenchmarkTools.BenchmarkGroup
@assert oneAPI_results isa BenchmarkTools.BenchmarkGroup

# Add other results to CPU results
for benchark in keys(CPU_results)
    for sim_method in keys(CPU_results[benchmark])
        CPU_results[benchmark][sim_method]["GPU"]["AMDGPU"] = AMD_results[benchmark][sim_method]["GPU"]["AMDGPU"]
        CPU_results[benchmark][sim_method]["GPU"]["CUDA"] = CUDA_results[benchmark][sim_method]["GPU"]["CUDA"]
        CPU_results[benchmark][sim_method]["GPU"]["Metal"] = Metal_results[benchmark][sim_method]["GPU"]["Metal"]
        CPU_results[benchmark][sim_method]["GPU"]["oneAPI"] = oneAPI_results[benchmark][sim_method]["GPU"]["oneAPI"]
    end
end

BenchmarkTools.save(joinpath(dirname(@__FILE__), "results", "combinedbenchmarks.json"), CPU_results)
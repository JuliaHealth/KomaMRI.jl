using BenchmarkTools

AMDGPU_results = nothing 
CUDA_results = nothing
Metal_results = nothing
oneAPI_results = nothing
const CPU_filepath = joinpath(dirname(@__FILE__), "results", "CPUbenchmarks.json")
const AMDGPU_filepath = joinpath(dirname(@__FILE__), "results", "AMDGPUbenchmarks.json")
const CUDA_filepath = joinpath(dirname(@__FILE__), "results", "CUDAbenchmarks.json")
const Metal_filepath = joinpath(dirname(@__FILE__), "results", "Metalbenchmarks.json")
const oneAPI_filepath = joinpath(dirname(@__FILE__), "results", "oneAPIbenchmarks.json")

if ispath(CPU_filepath) CPU_results = BenchmarkTools.load(CPU_filepath)[1] end
if ispath(AMDGPU_filepath) AMDGPU_results = BenchmarkTools.load(AMDGPU_filepath)[1] end
if ispath(CUDA_filepath) CUDA_results = BenchmarkTools.load(CUDA_filepath)[1] end
if ispath(Metal_filepath) Metal_results = BenchmarkTools.load(Metal_filepath)[1] end
if ispath(oneAPI_filepath) oneAPI_results = BenchmarkTools.load(oneAPI_filepath)[1] end

# Add other results to CPU results
@assert CPU_results isa BenchmarkTools.BenchmarkGroup
for benchmark in keys(CPU_results)
    for sim_method in keys(CPU_results[benchmark])
        if AMDGPU_results isa BenchmarkTools.BenchmarkGroup CPU_results[benchmark][sim_method]["GPU"]["AMDGPU"] = AMDGPU_results[benchmark][sim_method]["GPU"]["AMDGPU"] end
        if CUDA_results isa BenchmarkTools.BenchmarkGroup CPU_results[benchmark][sim_method]["GPU"]["CUDA"] = CUDA_results[benchmark][sim_method]["GPU"]["CUDA"] end
        if Metal_results isa BenchmarkTools.BenchmarkGroup CPU_results[benchmark][sim_method]["GPU"]["Metal"] = Metal_results[benchmark][sim_method]["GPU"]["Metal"] end
        if oneAPI_results isa BenchmarkTools.BenchmarkGroup CPU_results[benchmark][sim_method]["GPU"]["oneAPI"] = oneAPI_results[benchmark][sim_method]["GPU"]["oneAPI"] end
    end
end

BenchmarkTools.save(joinpath(dirname(@__FILE__), "results", "combinedbenchmarks.json"), CPU_results)
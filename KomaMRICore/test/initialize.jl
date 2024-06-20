using Pkg
using Suppressor

const USE_GPU = if "AMDGPU" in ARGS
    @suppress Pkg.add("AMDGPU")
    using AMDGPU
    @info "Testing AMD" maxlog=1
    true
elseif "CUDA" in ARGS
    @suppress Pkg.add("CUDA")
    using CUDA
    @info "Testing CUDA" maxlog=1
    true
elseif "Metal" in ARGS
    @suppress Pkg.add("Metal")
    using Metal
    @info "Testing Metal" maxlog=1
    true
elseif "oneAPI" in ARGS
    @suppress Pkg.add("oneAPI")
    using oneAPI
    @info "Testing oneAPI" maxlog=1
    true
else
    @info "Testing on the CPU with $(Threads.nthreads()) thread(s)" maxlog=1
    false
end

using KomaMRICore
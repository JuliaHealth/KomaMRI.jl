using Pkg
using Suppressor

const USE_GPU = if "AMDGPU" in ARGS
    using AMDGPU # ] add AMDGPU to KomaMRICore/test/Project.toml 
    @info "Testing AMD" maxlog=1
    true
elseif "CUDA" in ARGS
    using CUDA # ] add CUDA to KomaMRICore/test/Project.toml 
    @info "Testing CUDA" maxlog=1
    true
elseif "Metal" in ARGS
    using Metal # ] add Metal to KomaMRICore/test/Project.toml 
    @info "Testing Metal" maxlog=1
    true
elseif "oneAPI" in ARGS
    using oneAPI # ] add oneAPI to KomaMRICore/test/Project.toml 
    @info "Testing oneAPI" maxlog=1
    true
else
    @info "Testing on the CPU with $(Threads.nthreads()) thread(s)" maxlog=1
    false
end

using KomaMRICore
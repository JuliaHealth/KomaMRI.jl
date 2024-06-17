const USE_GPU = if "AMDGPU" in ARGS
    using AMDGPU
    @info "Testing AMD" maxlog=1
    true
elseif "CUDA" in ARGS
    using CUDA
    @info "Testing CUDA" maxlog=1
    true
elseif "Metal" in ARGS
    using Metal
    @info "Testing Metal" maxlog=1
    true
elseif "oneAPI" in ARGS
    using oneAPI
    @info "Testing oneAPI" maxlog=1
    true
else
    @info "Testing on the CPU with $(Threads.nthreads()) thread(s)" maxlog=1
    false
end

using KomaMRICore
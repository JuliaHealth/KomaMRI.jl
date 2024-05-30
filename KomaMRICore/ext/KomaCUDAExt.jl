module KomaCUDAExt

using CUDA
using KomaMRICore

name(::CUDABackend) = "CUDA"
isfunctional(::CUDABackend) = CUDA.functional()
set_device!(::CUDABackend, val) = CUDA.device!(val)
gpu_name(::CUDABackend) = CUDA.name(CUDA.device())
function reclaim_gpu(::CUDABackend)
    GC.gc(true)
    CUDA.reclaim()
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], CUDABackend())
end

end
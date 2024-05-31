module KomaCUDAExt

using CUDA
using KomaMRICore

name(::CUDABackend) = "CUDA"
isfunctional(::CUDABackend) = CUDA.functional()
set_device!(::CUDABackend, val) = CUDA.device!(val)
gpu_name(::CUDABackend) = CUDA.name(CUDA.device())

function print_gpus(::CUDABackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => CUDA.name(d) for
        (i, d) in enumerate(CUDA.devices())
    ]
    @info "$(length(CUDA.devices())) CUDA capable device(s)." devices...
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], CUDABackend())
end

end
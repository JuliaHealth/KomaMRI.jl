module KomaCUDAExt

using CUDA
import KomaMRICore
import Adapt

KomaMRICore.name(::CUDABackend) = "CUDA"
KomaMRICore.isfunctional(::CUDABackend) = CUDA.functional()
KomaMRICore.set_device!(::CUDABackend, val) = CUDA.device!(val)
KomaMRICore.device_name(::CUDABackend) = CUDA.name(CUDA.device())
@inline KomaMRICore._cis(x) = cis(x)

function KomaMRICore._print_devices(::CUDABackend)
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
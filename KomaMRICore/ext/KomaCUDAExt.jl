module KomaCUDAExt

using CUDA
import KomaMRICore

KomaMRICore.name(::CUDABackend) = "CUDA"
KomaMRICore.isfunctional(::CUDABackend) = CUDA.functional()
KomaMRICore.set_device!(::CUDABackend, val) = CUDA.device!(val)
KomaMRICore.device_name(::CUDABackend) = CUDA.name(CUDA.device())

function adapt_storage(::CUDABackend, x::ArbitraryMotion)
    fields = []
    for field in fieldnames(ArbitraryMotion)
        if field in (:ux, :uy, :uz) 
            push!(fields, adapt(::CUDABackend, getfield(x, field)))
        else
            push!(fields, f32(getfield(x, field)))
        end
    end
    return ArbitraryMotion(fields...)
end
function adapt_storage(
    ::CUDABackend, x::Vector{LinearInterpolator{T,V}}
) where {T<:Real,V<:AbstractVector{T}}
    return CUDA.cu.(x)
end

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
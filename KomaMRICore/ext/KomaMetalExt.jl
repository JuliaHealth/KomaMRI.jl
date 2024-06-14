module KomaMetalExt

using Metal
import KomaMRICore
import Adapt

KomaMRICore.name(::MetalBackend) = "Metal"
KomaMRICore.isfunctional(::MetalBackend) = Metal.functional()
KomaMRICore.set_device!(::MetalBackend, device_index::Integer) = device_index == 1 || @warn "Metal does not support multiple gpu devices. Ignoring the device setting."
KomaMRICore.set_device!(::MetalBackend, dev::Metal.MTLDevice) = Metal.device!(dev)
KomaMRICore.device_name(::MetalBackend) = String(Metal.current_device().name)

Adapt.adapt_storage(::MetalBackend, x::KomaMRICore.SimpleMotion) = KomaMRICore.f32(x)
function Adapt.adapt_storage(::MetalBackend, x::KomaMRICore.ArbitraryMotion)
    fields = []
    for field in fieldnames(KomaMRICore.ArbitraryMotion)
        if field in (:ux, :uy, :uz) 
            push!(fields, Adapt.adapt(MetalBackend(), getfield(x, field)))
        else
            push!(fields, KomaMRICore.f32(getfield(x, field)))
        end
    end
    return KomaMRICore.ArbitraryMotion(fields...)
end
function Adapt.adapt_storage(
    ::MetalBackend, x::Vector{KomaMRICore.LinearInterpolator{T,V}}
) where {T<:Real,V<:AbstractVector{T}}
    return Metal.mtl.(x)
end

function KomaMRICore._print_devices(::MetalBackend)
    @info "Metal device type: $(KomaMRICore.device_name(MetalBackend()))"
end

#Temporary workaround for https://github.com/JuliaGPU/Metal.jl/issues/348
#Once run_spin_excitation! and run_spin_precession! are kernel-based, this code
#can be removed
Base.cumsum(x::MtlVector) = convert(MtlVector, cumsum(KomaMRICore.cpu(x)))
Base.cumsum(x::MtlArray{T}; dims) where T = convert(MtlArray{T}, cumsum(KomaMRICore.cpu(x), dims=dims))
Base.findall(x::MtlVector{Bool}) = convert(MtlVector, findall(KomaMRICore.cpu(x)))

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], MetalBackend())
    @warn "Due to https://github.com/JuliaGPU/Metal.jl/issues/348, some functions may need to run on the CPU. Performance may be impacted as a result."
end

end
module KomaAMDGPUExt

using AMDGPU
import KomaMRICore
import Adapt

KomaMRICore.name(::ROCBackend) = "AMDGPU"
KomaMRICore.isfunctional(::ROCBackend) = AMDGPU.functional()
KomaMRICore.set_device!(::ROCBackend, dev_idx::Integer) = AMDGPU.device_id!(dev_idx)
KomaMRICore.set_device!(::ROCBackend, dev::AMDGPU.HIPDevice) = AMDGPU.device!(dev)
KomaMRICore.device_name(::ROCBackend) = AMDGPU.device().name

Adapt.adapt_storage(::ROCBackend, x::KomaMRICore.SimpleMotion) = KomaMRICore.f32(x)
function Adapt.adapt_storage(::ROCBackend, x::KomaMRICore.ArbitraryMotion)
    fields = []
    for field in fieldnames(KomaMRICore.ArbitraryMotion)
        if field in (:ux, :uy, :uz) 
            push!(fields, Adapt.adapt(ROCBackend(), getfield(x, field)))
        else
            push!(fields, KomaMRICore.f32(getfield(x, field)))
        end
    end
    return KomaMRICore.ArbitraryMotion(fields...)
end
function Adapt.adapt_storage(
    ::ROCBackend, x::Vector{KomaMRICore.LinearInterpolator{T,V}}
) where {T<:Real,V<:AbstractVector{T}}
    return AMDGPU.rocconvert.(x)
end

function KomaMRICore._print_devices(::ROCBackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => d.name for
        (i, d) in enumerate(AMDGPU.devices())
    ]
    @info "$(length(AMDGPU.devices())) AMD capable device(s)." devices...
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], ROCBackend())
end

end
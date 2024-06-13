module KomaAMDGPUExt

using AMDGPU
import KomaMRICore

KomaMRICore.name(::ROCBackend) = "AMDGPU"
KomaMRICore.isfunctional(::ROCBackend) = AMDGPU.functional()
KomaMRICore.set_device!(::ROCBackend, dev_idx::Integer) = AMDGPU.device_id!(dev_idx)
KomaMRICore.set_device!(::ROCBackend, dev::AMDGPU.HIPDevice) = AMDGPU.device!(dev)
KomaMRICore.device_name(::ROCBackend) = AMDGPU.device().name

function adapt_storage(::ROCBackend, x::ArbitraryMotion)
    fields = []
    for field in fieldnames(ArbitraryMotion)
        if field in (:ux, :uy, :uz) 
            push!(fields, adapt(::ROCBackend, getfield(x, field)))
        else
            push!(fields, f32(getfield(x, field)))
        end
    end
    return ArbitraryMotion(fields...)
end
function adapt_storage(
    ::ROCBackend, x::Vector{LinearInterpolator{T,V}}
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
module KomaoneAPIExt

using oneAPI
import KomaMRICore
import Adapt

KomaMRICore.name(::oneAPIBackend) = "oneAPI"
KomaMRICore.isfunctional(::oneAPIBackend) = oneAPI.functional()
KomaMRICore.set_device!(::oneAPIBackend, val) = oneAPI.device!(val)
KomaMRICore.device_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name

Adapt.adapt_storage(::oneAPIBackend, x::KomaMRICore.SimpleMotion) = KomaMRICore.f32(x)
function Adapt.adapt_storage(::oneAPIBackend, x::KomaMRICore.ArbitraryMotion)
    fields = []
    for field in fieldnames(KomaMRICore.ArbitraryMotion)
        if field in (:ux, :uy, :uz) 
            push!(fields, Adapt.adapt(oneAPIBackend(), getfield(x, field)))
        else
            push!(fields, KomaMRICore.f32(getfield(x, field)))
        end
    end
    return KomaMRICore.ArbitraryMotion(fields...)
end
function Adapt.adapt_storage(
    ::oneAPIBackend, x::Vector{KomaMRICore.LinearInterpolator{T,V}}
) where {T<:Real,V<:AbstractVector{T}}
    return oneAPI.oneArray.(x)
end

function KomaMRICore._print_devices(::oneAPIBackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => oneAPI.properties(d).name for
        (i, d) in enumerate(oneAPI.devices())
    ]
    @info "$(length(oneAPI.devices())) oneAPI capable device(s)." devices...
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], oneAPIBackend())
end

end
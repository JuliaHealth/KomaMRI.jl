module KomaoneAPIExt

using oneAPI
import KomaMRICore

KomaMRICore.name(::oneAPIBackend) = "oneAPI"
KomaMRICore.isfunctional(::oneAPIBackend) = oneAPI.functional()
KomaMRICore.set_device!(::oneAPIBackend, val) = oneAPI.device!(val)
KomaMRICore.device_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name

function adapt_storage(::oneAPIBackend, x::ArbitraryMotion)
    fields = []
    for field in fieldnames(ArbitraryMotion)
        if field in (:ux, :uy, :uz) 
            push!(fields, adapt(::oneAPIBackend, getfield(x, field)))
        else
            push!(fields, f32(getfield(x, field)))
        end
    end
    return ArbitraryMotion(fields...)
end
function adapt_storage(
    ::oneAPIBackend, x::Vector{LinearInterpolator{T,V}}
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
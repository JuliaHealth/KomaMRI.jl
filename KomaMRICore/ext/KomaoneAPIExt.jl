module KomaoneAPIExt

using oneAPI
import KomaMRICore
import Adapt

KomaMRICore.name(::oneAPIBackend) = "oneAPI"
KomaMRICore.isfunctional(::oneAPIBackend) = oneAPI.functional()
KomaMRICore.set_device!(::oneAPIBackend, val) = oneAPI.device!(val)
KomaMRICore.device_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name

function Adapt.adapt_storage(
    ::oneAPIBackend, x::Vector{KomaMRICore.LinearInterpolator{T,V}}
) where {T<:Real,V<:AbstractVector{T}}
    return Adapt.adapt.(oneArray, a)
end

function KomaMRICore._print_devices(::oneAPIBackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => oneAPI.properties(d).name for
        (i, d) in enumerate(oneAPI.devices())
    ]
    @info "$(length(oneAPI.devices())) oneAPI capable device(s)." devices...
end

#Temporary workaround since oneAPI.jl (similar to Metal) does not support some array operations
#Once run_spin_excitation! and run_spin_precession! are kernel-based, this code can be removed
Base.cumsum(x::oneVector) = convert(oneVector, cumsum(KomaMRICore.cpu(x)))
Base.cumsum(x::oneArray{T}; dims) where T = convert(oneArray{T}, cumsum(KomaMRICore.cpu(x), dims=dims))
Base.findall(x::oneVector{Bool}) = convert(oneVector, findall(KomaMRICore.cpu(x)))

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], oneAPIBackend())
    @warn "oneAPI does not support all array operations used by KomaMRI. GPU performance may be negatively impacted"
end

end
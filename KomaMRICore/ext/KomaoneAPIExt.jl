module KomaoneAPIExt

using oneAPI
import KomaMRICore
import Adapt

KomaMRICore.name(::oneAPIBackend) = "oneAPI"
KomaMRICore.isfunctional(::oneAPIBackend) = oneAPI.functional()
KomaMRICore.set_device!(::oneAPIBackend, val) = oneAPI.device!(val)
KomaMRICore.device_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name

KomaMRICore._devices(::oneAPIBackend) = [
    (id=i, name=oneAPI.properties(d).name, current=d == oneAPI.device()) for (i, d) in enumerate(oneAPI.devices())
]

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], oneAPIBackend())
    @warn "oneAPI support is experimental and does not support all array operations used by KomaMRI. GPU performance may be slower than expected"
end

end

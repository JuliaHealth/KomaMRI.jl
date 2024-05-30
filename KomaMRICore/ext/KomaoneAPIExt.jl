module KomaoneAPIExt

using oneAPI
using KomaMRICore

name(::oneAPIBackend) = "oneAPI"
isfunctional(::oneAPIBackend) = oneAPI.functional()
set_device(::oneAPIBackend, val) = oneAPI.device!(val)
gpu_name(::oneAPIBackend) = show(oneAPI.device())  

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], oneAPIBackend())
end

end
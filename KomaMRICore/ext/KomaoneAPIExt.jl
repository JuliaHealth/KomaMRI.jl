module KomaoneAPIExt

using oneAPI
using KomaMRICore

name(::oneAPIBackend) = "oneAPI"
isfunctional(::oneAPIBackend) = oneAPI.functional()
set_device(::oneAPIBackend, val) = oneAPI.device!(val)
gpu_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name

function print_gpus(::oneAPIBackend)
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
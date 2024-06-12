module KomaoneAPIExt

using oneAPI
import KomaMRICore

KomaMRICore.name(::oneAPIBackend) = "oneAPI"
KomaMRICore.isfunctional(::oneAPIBackend) = oneAPI.functional()
KomaMRICore.set_device!(::oneAPIBackend, val) = oneAPI.device!(val)
KomaMRICore.device_name(::oneAPIBackend) = oneAPI.properties(oneAPI.device()).name

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
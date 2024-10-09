## COV_EXCL_START

module KomaMetalExt

using Metal
import KomaMRICore
import Adapt

KomaMRICore.name(::MetalBackend) = "Metal"
KomaMRICore.isfunctional(::MetalBackend) = Metal.functional()
KomaMRICore.set_device!(::MetalBackend, device_index::Integer) = device_index == 1 || @warn "Metal does not support multiple gpu devices. Ignoring the device setting."
KomaMRICore.set_device!(::MetalBackend, dev::Metal.MTLDevice) = Metal.device!(dev)
KomaMRICore.device_name(::MetalBackend) = String(Metal.current_device().name)

function KomaMRICore._print_devices(::MetalBackend)
    @info "Metal device type: $(KomaMRICore.device_name(MetalBackend()))"
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], MetalBackend())
end

end

## COV_EXCL_STOP
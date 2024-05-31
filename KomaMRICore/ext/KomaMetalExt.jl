module KomaMetalExt

using Metal
using KomaMRICore

name(::MetalBackend) = "Metal"
isfunctional(::MetalBackend) = Metal.functional()
set_device!(::MetalBackend, device_index::Integer) = device_index == 1 || @warn "Metal does not support multiple gpu devices. Ignoring the device setting."
set_device!(::MetalBackend, dev::Metal.MTLDevice) = Metal.device!(dev)
gpu_name(::MetalBackend) = String(Metal.current_device().name)

function print_gpus(::MetalBackend)
    @info "Metal device type: $gpu_name(::MetalBackend)"
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], MetalBackend())
end

end
module KomaAMDGPUExt

using AMDGPU
using KomaMRICore

name(::ROCBackend) = "AMDGPU"
isfunctional(::ROCBackend) = AMDGPU.functional()
set_device(::ROCBackend, dev_idx::Integer) = AMDGPU.device_id!(dev_idx)
set_device(::ROCBackend, dev::AMDGPU.HIPDevice) = AMDGPU.device!(dev)
gpu_name(::ROCBackend) = AMDGPU.device().name

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], ROCBackend())
end

end
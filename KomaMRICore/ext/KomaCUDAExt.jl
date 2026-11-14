module KomaCUDAExt

using CUDA
using CUDA: NVML
import KomaMRICore
import Adapt

KomaMRICore.name(::CUDABackend) = "CUDA"
KomaMRICore.isfunctional(::CUDABackend) = CUDA.functional()
KomaMRICore.supports_warp_reduction(::CUDABackend) = true
KomaMRICore.set_device!(::CUDABackend, val) = CUDA.device!(val)
KomaMRICore.device_name(::CUDABackend) = KomaMRICore.device_name(CUDA.device())

function KomaMRICore.device_name(device::CUDA.CuDevice)
    CUDA.has_nvml() || return CUDA.name(device)
    id = NVML.index(NVML.Device(CUDA.parent_uuid(device)))
    return "$(CUDA.name(device)) — Device $id"
end

@inline function KomaMRICore.shfl_down(val, offset)
    CUDA.shfl_down_sync(0xffffffff, val, offset)
end

KomaMRICore._devices(::CUDABackend) = [
    (id=i-1, name=KomaMRICore.device_name(d), current=d == CUDA.device()) for (i, d) in enumerate(CUDA.devices())
]

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], CUDABackend())
end

"""Precompile CUDA GPU simulation workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

# CUDA.jl's persistent compiler cache requires Julia 1.11 or later.
@static if VERSION >= v"1.11"
    @setup_workload begin
        if CUDA.functional()
            KomaMRICore.BACKEND[] = CUDABackend()
            @compile_workload begin
                using KomaMRIBase
                using KomaMRICore

                sys = Scanner()
                obj_minimal = Phantom(
                    x=[0.0, 1e-3],
                    y=[0.0, 0.0],
                    z=[0.0, 0.0],
                    ρ=[1.0, 1.0],
                    T1=[0.8, 0.8],
                    T2=[80e-3, 80e-3],
                    T2s=[80e-3, 80e-3]
                )

                seq = PD.build_test_seq()

                sim_methods = [
                    Bloch(),
                    BlochMagnus2(),
                    BlochMagnus4(),
                ]
                precisions = ["f32"]
                return_types = ["mat", "raw"]

                for sim_method in sim_methods
                    for precision in precisions
                        for return_type in return_types
                            sim_params = Dict{String,Any}(
                                "sim_method" => sim_method,
                                "precision" => precision,
                                "return_type" => return_type,
                                "gpu" => true
                            )
                            signal = simulate(obj_minimal, seq, sys; sim_params, verbose=false)
                        end
                    end
                end
            end
            KomaMRICore.BACKEND[] = nothing
        end
    end
end

end

module KomaCUDAExt

using CUDA
import KomaMRICore
import Adapt

KomaMRICore.name(::CUDABackend) = "CUDA"
KomaMRICore.isfunctional(::CUDABackend) = CUDA.functional()
KomaMRICore.supports_warp_reduction(::CUDABackend) = true
KomaMRICore.set_device!(::CUDABackend, val) = CUDA.device!(val)
KomaMRICore.device_name(::CUDABackend) = CUDA.name(CUDA.device())

@inline function KomaMRICore.shfl_down(val, offset)
    CUDA.shfl_down_sync(0xffffffff, val, offset)
end

function KomaMRICore._print_devices(::CUDABackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => CUDA.name(d) for
        (i, d) in enumerate(CUDA.devices())
    ]
    @info "$(length(CUDA.devices())) CUDA capable device(s)." devices...
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], CUDABackend())
end

"""Precompile CUDA GPU simulation workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

# CUDA.jl's persistent compiler cache requires Julia 1.11 or later.
if VERSION >= v"1.11"
    @setup_workload begin
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

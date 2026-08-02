module KomaAMDGPUExt

using AMDGPU
import KomaMRICore
import Adapt

KomaMRICore.name(::ROCBackend) = "AMDGPU"
KomaMRICore.isfunctional(::ROCBackend) = AMDGPU.functional()
KomaMRICore.set_device!(::ROCBackend, dev_idx::Integer) = AMDGPU.device_id!(dev_idx)
KomaMRICore.set_device!(::ROCBackend, dev::AMDGPU.HIPDevice) = AMDGPU.device!(dev)
KomaMRICore.device_name(::ROCBackend) = AMDGPU.HIP.name(AMDGPU.device())

function KomaMRICore._print_devices(::ROCBackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => AMDGPU.HIP.name(d) for
        (i, d) in enumerate(AMDGPU.devices())
    ]
    @info "$(length(AMDGPU.devices())) AMD capable device(s)." devices...
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], ROCBackend())
end

"""Precompile AMDGPU simulation workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
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
        
        sim_methods = [Bloch(), BlochSimple()]
        precisions = ["f32", "f64"]
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
end

end

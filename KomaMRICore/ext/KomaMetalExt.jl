## COV_EXCL_START

module KomaMetalExt

using Metal
import KomaMRICore
import Adapt

KomaMRICore.name(::MetalBackend) = "Metal"
KomaMRICore.isfunctional(::MetalBackend) = Metal.functional()
KomaMRICore.supports_warp_reduction(::MetalBackend) = true
KomaMRICore.set_device!(::MetalBackend, device_index::Integer) = device_index == 1 || @warn "Metal does not support multiple gpu devices. Ignoring the device setting."
KomaMRICore.set_device!(::MetalBackend, dev::Metal.MTLDevice) = Metal.device!(dev)
KomaMRICore.device_name(::MetalBackend) = String(Metal.device().name)

@inline function KomaMRICore.shfl_down(val, offset)
    Metal.simd_shuffle_down(val, offset) 
end

function KomaMRICore._print_devices(::MetalBackend)
    @info "Metal device type: $(KomaMRICore.device_name(MetalBackend()))"
end

function __init__()
    push!(KomaMRICore.LOADED_BACKENDS[], MetalBackend())
end

end


"""Precompile Metal GPU simulation workflows for reduced first-use latency."""

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
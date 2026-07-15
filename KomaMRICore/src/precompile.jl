"""Precompile common KomaMRICore simulation workflows for reduced first-use latency."""

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
        
        seq = Sequence()
        append!(seq, PD.build_trapezoid(:x; area=8e-6, sys))
        append!(seq, PD.build_block_pulse(pi/2; duration=1e-3, sys, use=Excitation()))
        append!(seq, PD.build_sinc_pulse(pi/2; duration=2e-3, sys, use=Excitation()))
        append!(seq, PD.build_arbitrary_rf([1, 2, 1], pi/4; dwell=200e-6, sys, use=Excitation()))
        append!(seq, PD.build_adc(16; dwell=1e-6, sys))
        append!(seq, PD.build_label(:SET, :LIN, 3; sys))
        append!(seq, PD.build_rotation(pi/3; sys))
        append!(seq, PD.build_trigger(:physio1; sys))
        append!(seq, PD.build_arbitrary_grad(:x, [0, 1, 0]; sys))
        
        sim_methods = [
            Bloch(),
            BlochSimple(),
            BlochMagnusConst1(),
            BlochMagnusLin2(),
            BlochMagnusMid2(),
            BlochMagnusQuad2(),
            BlochMagnusQuad4(),
            BlochMagnusGL4(),
        ]
        
        precisions = ["f32", "f64"]
        return_types = ["mat", "raw"]
        
        for sim_method in sim_methods
            for precision in precisions
                for return_type in return_types
                    sim_params = Dict{String,Any}(
                        "sim_method" => sim_method,
                        "precision" => precision,
                        "return_type" => return_type,
                        "gpu" => false
                    )
                    signal = simulate(obj_minimal, seq, sys; sim_params, verbose=false)
                end
            end
        end
    end
end
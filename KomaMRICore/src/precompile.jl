"""Precompile common KomaMRICore simulation workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
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
            append!(seq, PD.build_trapezoid(:x; area=8e-6u"T*s/m", sys))
            append!(seq, PD.build_block_pulse(90u"deg"; duration=1u"ms", sys, use=Excitation()))
            append!(seq, PD.build_sinc_pulse(90u"deg"; duration=2u"ms", sys, use=Excitation()))
            append!(seq, PD.build_arbitrary_rf([1, 2, 1], 45u"deg"; dwell=200u"μs", sys, use=Excitation()))
            append!(seq, PD.build_adc(16; dwell=1u"μs", sys))
            append!(seq, PD.build_label(:SET, :LIN, 3; sys))
            append!(seq, PD.build_rotation(60u"deg"; sys))
            append!(seq, PD.build_trigger(:physio1; sys))
            append!(seq, PD.build_arbitrary_grad(:x, [0, 1, 0] .* u"mT/m"; sys))
            
            sim_methods = [
                Bloch(),
                BlochSimple(),
                BlochDict(),
                BlochMagnusConst1(),
                BlochMagnusLin2(),
                BlochMagnusMid2(),
                BlochMagnusLinComm2(),
                BlochMagnusQuad2(),
                BlochMagnusQuad4(),
                BlochMagnusGL2(),
                BlochMagnusGL4(),
                BlochMagnusBGL4(),
                BlochMagnusBGL6(),
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
                        signal = simulate(obj_minimal, seq, sys; sim_params)
                    end
                end
            end
        end
    end
end
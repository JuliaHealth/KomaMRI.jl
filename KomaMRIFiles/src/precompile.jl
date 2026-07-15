"""Precompile common KomaMRIFiles I/O workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        using KomaMRIBase
        using KomaMRICore
        
        sys = Scanner()
        obj = brain_phantom2D()
        
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
        
        mktempdir() do tmpdir
            seq_path = joinpath(tmpdir, "test.seq")
            write_seq(seq, seq_path; sys=sys)
            seq_read = read_seq(seq_path)
            
            phantom_path = joinpath(tmpdir, "phantom.phantom")
            write_phantom(obj, phantom_path)
            obj_read = read_phantom(phantom_path)
        end
    end
end
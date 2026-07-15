"""Precompile common KomaMRIFiles I/O workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
            using KomaMRIBase
            using KomaMRICore
            using Tempdir
            
            sys = Scanner()
            obj = brain_phantom2D()
            
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
end
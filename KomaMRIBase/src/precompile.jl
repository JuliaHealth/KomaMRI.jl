"""Precompile common KomaMRIBase workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
            sys = Scanner()
            obj = brain_phantom2D()
            
            seq = Sequence()
            
            PD.build_trapezoid(:x; area=8e-6u"T*s/m", sys)
            PD.build_block_pulse(90u"deg"; duration=1u"ms", sys, use=Excitation())
            PD.build_sinc_pulse(90u"deg"; duration=2u"ms", sys, use=Excitation())
            PD.build_arbitrary_rf([1, 2, 1], 45u"deg"; dwell=200u"μs", sys, use=Excitation())
            PD.build_adc(16; dwell=1u"μs", sys)
            PD.build_label(:SET, :LIN, 3; sys)
            PD.build_rotation(60u"deg"; sys)
            PD.build_trigger(:physio1; sys)
            PD.build_arbitrary_grad(:x, [0, 1, 0] .* u"mT/m"; sys)
            
            seq2 = Sequence()
            append!(seq2, PD.build_trapezoid(:x; area=8e-6u"T*s/m", sys))
            append!(seq2, PD.build_block_pulse(90u"deg"; duration=1u"ms", sys, use=Excitation()))
            append!(seq2, PD.build_sinc_pulse(90u"deg"; duration=2u"ms", sys, use=Excitation()))
            append!(seq2, PD.build_arbitrary_rf([1, 2, 1], 45u"deg"; dwell=200u"μs", sys, use=Excitation()))
            append!(seq2, PD.build_adc(16; dwell=1u"μs", sys))
            append!(seq2, PD.build_label(:SET, :LIN, 3; sys))
            append!(seq2, PD.build_rotation(60u"deg"; sys))
            append!(seq2, PD.build_trigger(:physio1; sys))
            append!(seq2, PD.build_arbitrary_grad(:x, [0, 1, 0] .* u"mT/m"; sys))
            
            seq_disc = discretize(seq2, sys)
        end
    end
end
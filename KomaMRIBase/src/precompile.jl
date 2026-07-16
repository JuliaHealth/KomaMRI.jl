"""Precompile common KomaMRIBase workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        sys = Scanner()
        obj = brain_phantom2D()
        
        seq = Sequence()
        @addblock begin
            seq += PD.build_trapezoid(:x; area=8e-6, sys)
            seq += PD.build_block_pulse(pi/2; duration=1e-3, sys, use=Excitation())
            seq += PD.build_sinc_pulse(pi/2; duration=2e-3, sys, use=Excitation())
            seq += PD.build_arbitrary_rf([1, 2, 1], pi/4; dwell=200e-6, sys, use=Excitation())
            seq += PD.build_adc(16; dwell=1e-6, sys)
            seq += PD.build_label(:SET, :LIN, 3; sys)
            seq += PD.build_rotation(pi/3; sys)
            seq += PD.build_trigger(:physio1; sys)
            seq += PD.build_arbitrary_grad(:x, [0, 1, 0]; sys)
        end
        
        seq_disc = discretize(seq, sys)
    end
end
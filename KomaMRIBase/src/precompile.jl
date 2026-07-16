"""Precompile common KomaMRIBase workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        sys = Scanner()
        obj = brain_phantom2D()
        
        seq = PD.build_test_seq()
        seq_disc = discretize(seq, sys)
    end
end
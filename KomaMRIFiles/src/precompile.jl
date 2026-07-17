"""Precompile common KomaMRIFiles I/O workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload
import KomaMRIBase: PulseDesigner as PD

@setup_workload begin
    @compile_workload begin
        using KomaMRIBase
        
        sys = Scanner()
        obj = brain_phantom2D()
        
        seq = PD.build_test_seq()
        
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

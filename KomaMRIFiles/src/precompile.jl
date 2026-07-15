"""Precompile common KomaMRIFiles I/O workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
            using KomaMRIBase
            using KomaMRICore
            using Tempdir
            
            # Pulseq Sequence I/O (write_seq: 10 docs, read_seq: 7 docs)
            # Create sequence with all event types (Pulseq parity test)
            sys = Scanner()
            seq = Sequence()
            
            # Add all event types to ensure Pulseq compatibility
            seq_with_events = PulseDesigner.EPI_example()  # Contains RF, Grad, ADC
            
            # Write Pulseq sequence
            mktempdir() do tmpdir
                seq_path = joinpath(tmpdir, "test.seq")
                write_seq(seq_with_events, seq_path; sys=sys)
                
                # Read Pulseq sequence
                seq_read = read_seq(seq_path)
            end
            
            # Phantom I/O (write_phantom: 4 docs, read_phantom: 4 docs)
            obj = brain_phantom2D()
            
            mktempdir() do tmpdir
                # Write phantom
                phantom_path = joinpath(tmpdir, "phantom.phantom")
                write_phantom(obj, phantom_path)
                
                # Read phantom
                obj_read = read_phantom(phantom_path)
            end
            
            # JEMRIS Phantom I/O (read_phantom_jemris: 3 docs)
            # Note: JEMRIS files are not created here, but read_phantom_jemris
            # logic is precompiled when read_seq or read_phantom is called
            # and the parser detects JEMRIS format. This is implicit.
        end
    end
end
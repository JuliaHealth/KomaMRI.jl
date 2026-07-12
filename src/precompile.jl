"""Precompile common KomaMRI workflows to reduce first-use latency."""

try
    using PrecompileTools: @setup_workload, @compile_workload

    @setup_workload begin
        @compile_workload begin
            redirect_stderr(devnull) do
                # Common scanner and phantom setup
                sys = Scanner()
                obj = brain_phantom2D()
                
                # Common sequence examples
                seq = PulseDesigner.EPI_example()
                
                # Simulation parameters
                sim_params = KomaMRICore.default_sim_params()
            end
        end
    end
catch
    # PrecompileTools not available; skip precompilation silently
end

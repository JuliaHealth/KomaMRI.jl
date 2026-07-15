"""Precompile common KomaMRIBase workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
            # Scanner Construction (38 docs, high frequency)
            sys = Scanner()
            
            # Phantom Construction (110 docs, highest frequency)
            # Built-in phantom
            obj = brain_phantom2D()
            
            # Custom phantom construction (different tissue types)
            x = collect(0.0:10e-3:20e-3)
            y = collect(0.0:10e-3:20e-3)
            z = collect(0.0:10e-3:20e-3)
            rho = ones(length(x))
            T1 = ones(length(x)) .* 0.8
            T2 = ones(length(x)) .* 80e-3
            T2s = T2
            obj_custom = Phantom(; x, y, z, ρ=rho, T1, T2, T2s)
            
            # Phantom with motion (NoMotion is most common)
            obj.motion = NoMotion()
            
            # Sequence Construction (178 docs, highest frequency)
            # Built-in sequence example
            seq = PulseDesigner.EPI_example()
            
            # Basic sequence with all event types
            seq_custom = Sequence()
            
            # Discretize sequence (critical for simulation)
            # Note: discretize may be in KomaMRICore, but include here
            # if it's part of base sequence operations
            
            # Sequence Event Types (all required for Pulseq parity)
            # RF event
            rf = RF(; A=1.0, T=1e-3)
            
            # Gradient events (x, y, z)
            gx = Grad(; A=1.0, T=1e-3)
            gy = Grad(; A=1.0, T=1e-3)
            gz = Grad(; A=1.0, T=1e-3)
            
            # ADC event
            adc = ADC(; N=100, T=10e-3, delay=0.0)
            
            # Delay event
            del = Delay(; T=1e-3)
        end
    end
end
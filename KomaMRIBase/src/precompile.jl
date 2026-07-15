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
            
            # RF Event Types (all variants required for Pulseq parity)
            # Block pulse (A and T are numbers)
            rf_block = RF(10e-3, 0.5e-3, 0, 0.1e-3)
            
            # Uniformly-sampled waveform (A is vector, T is number)
            tl = -3:0.2:-0.2
            tr = 0.2:0.2:3
            A_uniform = (10e-3) * [sin.(π*tl)./(π*tl); 1; sin.(π*tr)./(π*tr)]
            rf_uniform = RF(A_uniform, 0.5e-3, 0, 0.1e-3)
            
            # Time-shaped waveform (A and T are both vectors)
            tl = -4:0.2:-0.2
            tr = 0.2:0.2:4
            A_shaped = (10e-3) * [sin.(π*tl)./(π*tl); 1; 1; sin.(π*tr)./(π*tr)]
            T_shaped = [0.05e-3*ones(length(tl)); 2e-3; 0.05e-3*ones(length(tl))]
            rf_shaped = RF(A_shaped, T_shaped, 0, 0.1e-3)
            
            # Frequency modulated RF (with Δf)
            rf_fm = RF(10e-3, 0.5e-3, 1000, 0.1e-3)  # 1kHz offset
            
            # Gradient Event Types (all variants for x, y, z axes)
            # Trapezoidal gradient (A and T are numbers)
            gr_trap = Grad(50e-6, 5e-3, 1e-3, 1e-3, 2e-3)
            
            # Uniformly-sampled gradient (A is vector, T is number)
            t_grad = 0:0.25:7.5
            A_grad_uniform = 10e-6 * sqrt.(π*t_grad) .* sin.(π*t_grad)
            gr_uniform = Grad(A_grad_uniform, 10e-3, 0, 1e-3, 1e-3)
            
            # Time-shaped gradient (A and T are both vectors)
            A_grad_shaped = 50e-6 * [1; 1; 0.8; 0.8; 1; 1]
            T_grad_shaped = 1e-3 * [5; 0.2; 5; 0.2; 5]
            gr_shaped = Grad(A_grad_shaped, T_grad_shaped, 1e-3, 1e-3, 1e-3)
            
            # Add all gradient variants to sequence (x, y, z axes)
            seq_grad = Sequence()
            @addblock seq_grad += (x=gr_trap, y=gr_uniform, z=gr_shaped)
            
            # ADC Event
            adc = ADC(16, 5e-3, 1e-3)
            
            # Extension Events (Labels, Triggers, Rotations)
            # Label increment
            lInc = LabelInc(1, "LIN")
            
            # Label set
            lSet = LabelSet(1, "ECO")
            
            # Trigger extension
            trig = Trigger(0, 1, 100e-6, 500e-6)
            
            # Quaternion rotation extension
            qrot = QuaternionRot(1.0, 0.0, 0.0, 0.0)
            
            # Complete Sequence with All Event Types
            seq_complete = Sequence()
            @addblock seq_complete += (rf_block, adc, x=gr_trap, y=gr_uniform, z=gr_shaped, lInc, trig)
            @addblock seq_complete += (rf_uniform, x=gr_trap, lSet)
            @addblock seq_complete += (rf_shaped, adc, x=gr_uniform, y=gr_shaped)
        end
    end
end
"""Precompile common KomaMRICore simulation workflows for reduced first-use latency."""

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        redirect_stderr(devnull) do
            # Load base types
            using KomaMRIBase
            
            # Default Simulation Parameters (7 docs)
            sim_params_default = default_sim_params()
            
            # Simulate with All Simulation Methods (CRITICAL - 22 docs)
            # Use minimal 1-2 spin phantom for speed
            sys = Scanner()
            obj_minimal = Phantom(
                x=[0.0, 1e-3],
                y=[0.0, 0.0],
                z=[0.0, 0.0],
                ρ=[1.0, 1.0],
                T1=[0.8, 0.8],
                T2=[80e-3, 80e-3],
                T2s=[80e-3, 80e-3]
            )
            seq_minimal = PulseDesigner.EPI_example()
            
            # Bloch simulation method (most common)
            sim_params_bloch = default_sim_params()
            sim_params_bloch["sim_method"] = Bloch()
            signal_bloch = simulate(obj_minimal, seq_minimal, sys; sim_params=sim_params_bloch)
            
            # BlochDict simulation method (GPU-friendly)
            sim_params_dict = default_sim_params()
            sim_params_dict["sim_method"] = BlochDict()
            signal_dict = simulate(obj_minimal, seq_minimal, sys; sim_params=sim_params_dict)
            
            # Data Type Variants (Float32 and Float64)
            # Float32 (default, common for performance)
            sim_params_f32 = default_sim_params()
            sim_params_f32["return_type"] = "mat"  # Dense matrix
            signal_f32 = simulate(obj_minimal, seq_minimal, sys; sim_params=sim_params_f32)
            
            # Float64 (precision-critical applications)
            obj_f64 = Phantom(
                x=Float64[0.0, 1e-3],
                y=Float64[0.0, 0.0],
                z=Float64[0.0, 0.0],
                ρ=Float64[1.0, 1.0],
                T1=Float64[0.8, 0.8],
                T2=Float64[80e-3, 80e-3],
                T2s=Float64[80e-3, 80e-3]
            )
            signal_f64 = simulate(obj_f64, seq_minimal, sys; sim_params=sim_params_f32)
            
            # Export Formats (mat and raw)
            # Mat export (MATLAB-compatible)
            sim_params_mat = default_sim_params()
            sim_params_mat["return_type"] = "mat"
            signal_mat = simulate(obj_minimal, seq_minimal, sys; sim_params=sim_params_mat)
            
            # Raw export (native Julia types)
            sim_params_raw = default_sim_params()
            sim_params_raw["return_type"] = "raw"
            signal_raw = simulate(obj_minimal, seq_minimal, sys; sim_params=sim_params_raw)
            
            # GPU Path (if available)
            # This will be a no-op on CPU, but precompiles GPU code if available
            try
                sim_params_gpu = default_sim_params()
                sim_params_gpu["gpu"] = true
                signal_gpu = simulate(obj_minimal, seq_minimal, sys; sim_params=sim_params_gpu)
            catch
                # GPU not available; skip silently
            end
            
            # Discretize Sequence (critical for simulation prep)
            # Discretization with different methods
            seq_disc = discretize(seq_minimal, sys)
        end
    end
end
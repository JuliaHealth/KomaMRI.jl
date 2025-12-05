function setup_3DSlice_benchmark()
    sys = Scanner()
    obj = brain_phantom3D()[1:10000]
    seq = PulseDesigner.EPI_example()
    return (sys, obj, seq)
end

function setup_MRILab_benchmark()
    sys = Scanner()
    sys.Smax = 150    # [mT/m/ms]
    sys.Gmax = 500e-3 # [T/m]
    sys.GR_Δt = 4e-6  # [s]
    FOV = 0.2         # [m]
    N = 80            # Reconstructed image N×N
    ## Pulse programming
    # RF sinc
    B1 = 1e-6 # For 90 deg flip angle
    Trf = 1e-3
    rf = PulseDesigner.RF_sinc(B1, Trf, sys; TBP=4, a=0.46)
    α_desired = 90 + 0im
    α =  get_flip_angles(rf)[1]
    rf *= α_desired / α #Linearly adjusts B1 to obtain desired FA
    # Spiral sequence
    TE = 50e-3  # 50e-3 [s]
    TR = 10     # 10 [s]
    Nint = 8
    spiral = PulseDesigner.spiral_base(FOV, N, sys; BW=60e3, Nint, λ=2.1)
    delayTE = Delay(TE - Trf / 2)
    delayTR = Delay(TR - Trf / 2 - TE - dur(spiral(0)))
    seq = Sequence()
    for i = 1:Nint
        seq += rf + delayTE + spiral(i - 1) + delayTR
    end
    obj = brain_phantom2D()

    return (sys, obj, seq)
end

function add_benchmark_to_suite(suite::BenchmarkGroup, benchmark_name::String, backend::String, num_cpu_threads::Int64, sys::Scanner, obj::Phantom, seq::Sequence, sim_params::Dict)
    if backend == "CPU"
        suite[benchmark_name]["Bloch"]["CPU"][string(num_cpu_threads, " ", "thread(s)")] = @benchmarkable @suppress simulate($obj, $seq, $sys; sim_params=$sim_params)
    else
        suite[benchmark_name]["Bloch"]["GPU"][backend] = @benchmarkable @suppress simulate($obj, $seq, $sys; sim_params=$sim_params)
    end
end

function setup_benchmarks(suite::BenchmarkGroup, backend::String, num_cpu_threads::Int64)    
    # Benchmark 1: from lit-04-3DSliceSelective.jl
    sys1, obj1, seq1 = setup_3DSlice_benchmark()
    sim_params1 = Dict{String,Any}(
        "gpu" => (backend != "CPU"),
        "Nthreads" => num_cpu_threads
    )
    add_benchmark_to_suite(suite, "Slice Selection 3D", backend, num_cpu_threads, sys1, obj1, seq1, sim_params1)

    # Benchmark 2: from MRiLab_speed.jl
    sys2, obj2, seq2 = setup_MRILab_benchmark()
    sim_params2 = Dict{String,Any}(
        "gpu" => (backend != "CPU"),
        "Nthreads" => num_cpu_threads
    )
    add_benchmark_to_suite(suite, "MRI Lab", backend, num_cpu_threads, sys2, obj2, seq2, sim_params2)
end
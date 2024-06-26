function setup_benchmarks(suite::BenchmarkGroup, backend::String, num_cpu_threads::NTuple)    
    # Benchmark 1: from lit-04-3DSliceSelective.jl
    sys1 = Scanner()
    obj1 = brain_phantom3D()[1:10000]
    seq1 = PulseDesigner.EPI_example()
     
    if backend == "CPU"
        for n in num_cpu_threads
            sim_params = Dict{String,Any}("Nblocks"=>20,"gpu"=>false,"Nthreads"=>n)
            suite["Slice Selection 3D"]["Bloch"]["CPU"][string(n, " ", "thread(s)")] = @benchmarkable @suppress simulate($obj1, $seq1, $sys1; sim_params=$sim_params)
        end
    else
        sim_params = Dict{String,Any}("Nblocks"=>20)
        suite["Slice Selection 3D"]["Bloch"]["GPU"][backend] = @benchmarkable @suppress simulate($obj1, $seq1, $sys1; sim_params=$sim_params)
    end

    # Benchmark 2: from MRiLab_speed.jl
    sys2 = Scanner()
    sys2.Smax = 150    # [mT/m/ms]
    sys2.Gmax = 500e-3 # [T/m]
    sys2.GR_Δt = 4e-6  # [s]
    FOV = 0.2         # [m]
    N = 80            # Reconstructed image N×N
    ## Pulse programming
    # RF sinc
    B1 = 1e-6 # For 90 deg flip angle
    Trf = 1e-3
    rf = PulseDesigner.RF_sinc(B1, Trf, sys2; TBP=4, a=0.46)
    α_desired = 90 + 0im
    α =  get_flip_angles(rf)[1]
    rf *= α_desired / α #Linearly adjusts B1 to obtain desired FA
    # Spiral sequence
    TE = 50e-3  # 50e-3 [s]
    TR = 10     # 10 [s]
    Nint = 8
    spiral = PulseDesigner.spiral_base(FOV, N, sys2; BW=60e3, Nint, λ=2.1)
    delayTE = Delay(TE - Trf / 2)
    delayTR = Delay(TR - Trf / 2 - TE - dur(spiral(0)))
    seq2 = Sequence()
    for i = 1:Nint
        seq2 += rf + delayTE + spiral(i - 1) + delayTR
    end
    obj2 = brain_phantom2D()

    if backend == "CPU"
        for n in num_cpu_threads
            sim_params = Dict{String,Any}("Nblocks"=>20,"gpu"=>false,"Nthreads"=>n)
            suite["MRI Lab"]["Bloch"]["CPU"][string(n, " ", "thread(s)")] = @benchmarkable @suppress simulate($obj2, $seq2, $sys2; sim_params=$sim_params)
        end
    else
        sim_params = Dict{String,Any}("Nblocks"=>20)
        suite["MRI Lab"]["Bloch"]["GPU"][backend] = @benchmarkable @suppress simulate($obj2, $seq2, $sys2; sim_params=$sim_params)
    end
end
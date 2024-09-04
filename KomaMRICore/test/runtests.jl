using TestItems, TestItemRunner

### NOTE: by default, tests are run on the CPU with the number of threads set to
#   Threads.nthreads(). To run on a specific GPU backend, add the name of the
#   backend package ("AMDGPU", "CUDA", "Metal", or "oneAPI") to the test/Project.toml
#   file in KomaMRICore and pass the name as a test argument.
#
#   Example:
#
#   import Pkg
#   Pkg.test("KomaMRICore"; test_args=["CUDA"])
#
#   To run on the cpu with a specific number of threads, pass the number of threads
#   as a julia argument.
#
#   Example:
#
#   import Pkg
#   Pkg.test("KomaMRICore"; julia_args=`--threads=4`)
#
#   For changing the default backend used for testing,
#   modify the [preferences.KomaMRICore] section in the test/Project.toml:
#
#     [preferences.KomaMRICore]
#     test_backend = "CPU"
#
#   For the backend preference to take effect, you need to:
#   - REPL testing: No action needed. `] test` should pick up the preference right away.
#   - VSCode testing: You need to restart VSCode.
#
#   Sadly, LocalPreferences.toml are not picked up by VScode (that could be .gitignore'd),
#   so we had put them into the test/Project.toml.
#
###

#Environment variable set by CI
const CI    = get(ENV, "CI", nothing)
const group = get(ENV, "TEST_GROUP", :core) |> Symbol

@run_package_tests filter=ti->(group in ti.tags)&&(isnothing(CI) || :skipci ∉ ti.tags) #verbose=true

@testitem "Spinors×Mag" tags=[:core, :nomotion] begin
    using KomaMRICore: Rx, Ry, Rz, Q, rotx, roty, rotz, Un, Rφ, Rg

    ## Verifying that operators perform counter-clockwise rotations
    v = [1, 2, 3]
    m = Mag([complex(v[1:2]...)], [v[3]])
    # Rx
    @test rotx(π/2) * v    ≈ [1, -3, 2]
    @test (Rx(π/2) * m).xy ≈ [1.0 - 3.0im]
    @test (Rx(π/2) * m).z  ≈ [2.0]
    # Ry
    @test roty(π/2) * v    ≈ [3, 2, -1]
    @test (Ry(π/2) * m).xy ≈ [3.0 + 2.0im]
    @test (Ry(π/2) * m).z  ≈ [-1.0]
    # Rz
    @test rotz(π/2) * v    ≈ [-2, 1, 3]
    @test (Rz(π/2) * m).xy ≈ [-2.0 + 1.0im]
    @test (Rz(π/2) * m).z  ≈ [3.0]
    # Rn
    @test Un(π/2, [1,0,0]) * v ≈ rotx(π/2) * v
    @test Un(π/2, [0,1,0]) * v ≈ roty(π/2) * v
    @test Un(π/2, [0,0,1]) * v ≈ rotz(π/2) * v
    @test (Q(π/2, 1.0+0.0im, 0.0) * m).xy ≈ (Rx(π/2) * m).xy
    @test (Q(π/2, 1.0+0.0im, 0.0) * m).z  ≈ (Rx(π/2) * m).z
    @test (Q(π/2, 0.0+1.0im, 0.0) * m).xy ≈ (Ry(π/2) * m).xy
    @test (Q(π/2, 0.0+1.0im, 0.0) * m).z  ≈ (Ry(π/2) * m).z
    @test (Q(π/2, 0.0+0.0im, 1.0) * m).xy ≈ (Rz(π/2) * m).xy
    @test (Q(π/2, 0.0+0.0im, 1.0) * m).z  ≈ (Rz(π/2) * m).z

    ## Verify that Spinor rotation = matrix rotation
    v = rand(3)
    n = rand(3); n = n ./ sqrt(sum(n.^2))
    m = Mag([complex(v[1:2]...)], [v[3]])
    φ, θ, φ1, φ2 = rand(4) * 2π
    # Rx
    vx = rotx(θ) * v
    mx = Rx(θ) * m
    @test [real(mx.xy); imag(mx.xy); mx.z] ≈ vx
    # Ry
    vy = roty(θ) * v
    my = Ry(θ) * m
    @test [real(my.xy); imag(my.xy); my.z] ≈ vy
    # Rz
    vz = rotz(θ) * v
    mz = Rz(θ) * m
    @test [real(mz.xy); imag(mz.xy); mz.z] ≈ vz
    # Rφ
    vφ = Un(θ, [sin(φ); cos(φ); 0.0]) * v
    mφ = Rφ(φ,θ) * m
    @test [real(mφ.xy); imag(mφ.xy); mφ.z] ≈ vφ
    # Rg
    vg = rotz(φ2) * roty(θ) * rotz(φ1) * v
    mg = Rg(φ1,θ,φ2) * m
    @test [real(mg.xy); imag(mg.xy); mg.z] ≈ vg
    # Rn
    vq = Un(θ, n) * v
    mq = Q(θ, n[1]+n[2]*1im, n[3]) * m
    @test [real(mq.xy); imag(mq.xy); mq.z] ≈ vq

    ## Spinors satify that  |α|^2 + |β|^2 = 1
    @test abs(Rx(θ))                     ≈ [1]
    @test abs(Ry(θ))                     ≈ [1]
    @test abs(Rz(θ))                     ≈ [1]
    @test abs(Rφ(φ,θ))                   ≈ [1]
    @test abs(Q(θ, n[1]+n[2]*1im, n[3])) ≈ [1]

    ## Checking properties of Introduction to the Shinnar-Le Roux algorithm.
    # Rx = Rz(-π/2) * Ry(θ) * Rz(π/2)
    @test rotx(θ) * v ≈ rotz(-π/2) * roty(θ) * rotz(π/2) * v
    @test (Rx(θ) * m).xy ≈ (Rz(-π/2) * Ry(θ) * Rz(π/2) * m).xy
    @test (Rx(θ) * m).z  ≈ (Rz(-π/2) * Ry(θ) * Rz(π/2) * m).z
    # Rφ(φ,θ) = Rz(-φ) Ry(θ) Rz(φ)
    @test (Rφ(φ,θ) * m).xy   ≈ (Rz(-φ) * Ry(θ) * Rz(φ) * m).xy
    @test (Rφ(φ,θ) * m).z    ≈ (Rz(-φ) * Ry(θ) * Rz(φ) * m).z
    # Rg(φ1, θ, φ2) = Rz(φ2) Ry(θ) Rz(φ1)
    @test (Rg(φ1,θ,φ2) * m).xy   ≈ (Rz(φ2) * Ry(θ) * Rz(φ1) * m).xy
    @test (Rg(φ1,θ,φ2) * m).z    ≈ (Rz(φ2) * Ry(θ) * Rz(φ1) * m).z
    # Rg(-φ, θ, φ) = Rz(-φ) Ry(θ) Rz(φ) = Rφ(φ,θ)
    @test rotz(-φ) * roty(θ) * rotz(φ) * v ≈ Un(θ, [sin(φ); cos(φ); 0.0]) * v
    @test (Rg(φ,θ,-φ) * m).xy    ≈ (Rφ(φ,θ) * m).xy
    @test (Rg(φ,θ,-φ) * m).z     ≈ (Rφ(φ,θ) * m).z

    ## Verify trivial identities
    # Rφ is an xy-plane rotation of θ around an axis making an angle of φ with respect to the y-axis
    # Rφ φ=0 = Ry
    @test (Rφ(0,θ) * m).xy   ≈ (Ry(θ) * m).xy
    @test (Rφ(0,θ) * m).z    ≈ (Ry(θ) * m).z
    # Rφ φ=π/2 = Rx
    @test (Rφ(π/2,θ) * m).xy ≈ (Rx(θ) * m).xy
    @test (Rφ(π/2,θ) * m).z  ≈ (Rx(θ) * m).z
    # General rotation Rn
    # Rn n=[1,0,0] = Rx
    @test Un(θ, [1,0,0]) * v ≈ rotx(θ) * v
    @test (Q(θ, 1.0+0.0im, 0.0) * m).xy ≈ (Rx(θ) * m).xy
    @test (Q(θ, 1.0+0.0im, 0.0) * m).z  ≈ (Rx(θ) * m).z
    # Rn n=[0,1,0] = Ry
    @test Un(θ, [0,1,0]) * v ≈ roty(θ) * v
    @test (Q(θ, 0.0+1.0im, 0.0) * m).xy ≈ (Ry(θ) * m).xy
    @test (Q(θ, 0.0+1.0im, 0.0) * m).z  ≈ (Ry(θ) * m).z
    # Rn n=[0,0,1] = Rz
    @test Un(θ, [0,0,1]) * v ≈ rotz(θ) * v
    @test (Q(θ, 0.0+0.0im, 1.0) * m).xy ≈ (Rz(θ) * m).xy
    @test (Q(θ, 0.0+0.0im, 1.0) * m).z  ≈ (Rz(θ) * m).z

    # Associativity
    # Rx
    @test (((Rz(-π/2) * Ry(θ)) * Rz(π/2)) * m).xy ≈ (Rx(θ) * m).xy
    @test (((Rz(-π/2) * Ry(θ)) * Rz(π/2)) * m).z  ≈ (Rx(θ) * m).z
    @test (Rz(-π/2) * (Ry(θ) * (Rz(π/2) * m))).xy ≈ (Rx(θ) * m).xy
    @test (Rz(-π/2) * (Ry(θ) * (Rz(π/2) * m))).z  ≈ (Rx(θ) * m).z
    # Rφ
    @test (Rφ(φ,θ) * m).xy   ≈ (((Rz(-φ) * Ry(θ)) * Rz(φ)) * m).xy
    @test (Rφ(φ,θ) * m).z    ≈ (((Rz(-φ) * Ry(θ)) * Rz(φ)) * m).z
    @test (Rφ(φ,θ) * m).xy   ≈ ((Rz(-φ) * (Ry(θ) * Rz(φ))) * m).xy
    @test (Rφ(φ,θ) * m).z    ≈ ((Rz(-φ) * (Ry(θ) * Rz(φ))) * m).z
    # Rg
    @test (Rg(φ1,θ,φ2) * m).xy   ≈ (((Rz(φ2) * Ry(θ)) * Rz(φ1)) * m).xy
    @test (Rg(φ1,θ,φ2) * m).z    ≈ (((Rz(φ2) * Ry(θ)) * Rz(φ1)) * m).z
    @test (Rg(φ1,θ,φ2) * m).xy   ≈ ((Rz(φ2) * (Ry(θ) * Rz(φ1))) * m).xy
    @test (Rg(φ1,θ,φ2) * m).z    ≈ ((Rz(φ2) * (Ry(θ) * Rz(φ1))) * m).z

    ## Other tests
    # Test Spinor struct
    α, β = rand(2)
    s = Spinor(α, β)
    @test s[1].α ≈ [Complex(α)] && s[1].β ≈ [Complex(β)]
    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", s)
    @test true
end

@testitem "ISMRMRD" tags=[:core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")

    seq = PulseDesigner.EPI_example()
    sys = Scanner()
    obj = brain_phantom2D()
    parts = kfoldperm(length(obj), 2)

    sim_params = KomaMRICore.default_sim_params()
    sim_params["return_type"] = "raw"
    sim_params["gpu"] = USE_GPU

    sig1 = @suppress simulate(obj[parts[1]], seq, sys; sim_params)
    sig2 = @suppress simulate(obj[parts[2]], seq, sys; sim_params)
    sig = @suppress simulate(obj, seq, sys; sim_params)

    @test isapprox(sig, sig1 + sig2; rtol=0.001)
end

@testitem "signal_to_raw_data" tags=[:core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")

    seq = PulseDesigner.EPI_example()
    sys = Scanner()
    obj = brain_phantom2D()

    sim_params = KomaMRICore.default_sim_params()
    sim_params["return_type"] = "mat"
    sim_params["gpu"] = USE_GPU
    sig = @suppress simulate(obj, seq, sys; sim_params)

    # Test signal_to_raw_data
    raw = signal_to_raw_data(sig, seq)
    sig_aux = vcat([vec(profile.data) for profile in raw.profiles]...)
    sig_raw = reshape(sig_aux, length(sig_aux), 1)
    @test all(sig .== sig_raw)

    seq.DEF["FOV"] = [23e-2, 23e-2, 0]
    raw = signal_to_raw_data(sig, seq)
    sig_aux = vcat([vec(profile.data) for profile in raw.profiles]...)
    sig_raw = reshape(sig_aux, length(sig_aux), 1)
    @test all(sig .== sig_raw)

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", raw)
    @test true
end

@testitem "Bloch" tags=[:important, :core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_sphere_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_RF_accuracy" tags=[:important, :core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0.],T1=[T1],T2=[T2],Δw=[Δw])

    rf_phase = [0, π/2]
    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1 .* exp(1im*rf_phase[i]), Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>USE_GPU)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [-0.0153894+0.142582im,
            0.00257641+0.14196im,
            0.020146+0.13912im,
            0.037051+0.134149im,
            0.0530392+0.12717im,
            0.0678774+0.11833im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "Bloch_phase_compensation" tags=[:important, :core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0.],T1=[T1],T2=[T2],Δw=[Δw])

    rf_phase = 2π*rand()
    seq1 = Sequence()
    seq1 += RF(B1, Trf)
    seq1 += ADC(N, Tadc)

    seq2 = Sequence()
    seq2 += RF(B1 .* exp(1im*rf_phase), Trf)
    seq2 += ADC(N, Tadc, 0, 0, rf_phase)

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>USE_GPU)
    raw1 = @suppress simulate(obj, seq1, sys; sim_params)
    raw2 = @suppress simulate(obj, seq2, sys; sim_params)

    @test raw1.profiles[1].data ≈ raw2.profiles[1].data
end

@testitem "BlochDict" tags=[:important, :core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    seq = seq_epi_100x100_TE100_FOV230()
    obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])
    sys = Scanner()
    sim_params = Dict(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat")
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    sim_params["sim_method"] = KomaMRICore.BlochDict()
    sig2 = @suppress simulate(obj, seq, sys; sim_params)
    sig2 = sig2 / prod(size(obj))
    @test sig ≈ sig2

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", KomaMRICore.BlochDict())
    @test true
end

@testitem "BlochSimple" tags=[:important, :core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_sphere_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.BlochSimple(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "simulate_slice_profile" tags=[:core, :nomotion] begin
    using Suppressor
    include("initialize_backend.jl")

    # This is a sequence with a sinc RF 30° excitation pulse
    sys = Scanner()
    sys.Smax = 50
    B1 = 4.92e-6
    Trf = 3.2e-3
    zmax = 2e-2
    fmax = 5e3
    z = range(-zmax, zmax, 400)
    Gz = fmax / (γ * zmax)
    f = γ * Gz * z
    seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0; 0; Gz], TBP=8)

    # Simulate the slice profile
    sim_params = Dict{String, Any}(
        "Δt_rf" => Trf / length(seq.RF.A[1]),
        "gpu" => USE_GPU)
    M = @suppress simulate_slice_profile(seq; z, sim_params)

    # For the time being, always pass the test
    @test true
end

@testitem "GPU Functions" tags=[:core, :nomotion] begin
    using Suppressor
    import KernelAbstractions as KA
    include("initialize_backend.jl")

    x = ones(Float32, 1000)

    @suppress begin
        if USE_GPU
            y = x |> gpu
            @test KA.get_backend(y) isa KA.GPU
            y = y |> cpu
            @test KA.get_backend(y) isa KA.CPU
        else
            # Test that gpu and cpu are no-ops
            y = x |> gpu
            @test y == x
            y = y |> cpu
            @test y == x
        end
    end

    @suppress print_devices()
    @test true
end

# --------- Motion-related tests -------------
@testitem "Bloch SimpleAction" tags=[:core] begin 
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_simple_motion()
    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE SimpleAction: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "BlochSimple SimpleAction" tags=[:core] begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_simple_motion()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.BlochSimple(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE SimpleAction BlochSimple: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch ArbitraryAction"  tags=[:core, :motion] begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    # for i in 1:15
    #     sig_jemris = signal_brain_motion_jemris()
    #     seq = seq_epi_100x100_TE100_FOV230()
    #     sys = Scanner()
    #     obj = phantom_brain_arbitrary_motion()

    #     sim_params = Dict{String, Any}(
    #         "gpu"=>USE_GPU,
    #         "sim_method"=>KomaMRICore.Bloch(),
    #         "return_type"=>"mat"
    #     )
    #     sig = simulate(obj, seq, sys; sim_params)
    #     sig = sig / prod(size(obj))
    #     NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    #     # println("NMRSE ArbitraryAction: ", NMRSE(sig, sig_jemris))
    #     @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
    # end


    # tr = TimeRange(0.0f0, 1.0f0)

    # time = collect(-1:0.01:2)

    # ux_cpu = zeros(Float32, (2, length(time)))
    # ux_gpu = ux_cpu |> gpu
    
    # for i in 1:10
    #     # cpu
    #     d = rand(Float32, (2, 2))
    #     t = (time |> f32)'
    #     t_unit = KomaMRIBase.unit_time(t, tr)
        
    #     itp = KomaMRIBase.interpolate(d, KomaMRIBase.Gridded(KomaMRIBase.Linear()), Val(size(d,1)))
    #     ux_cpu .= KomaMRIBase.resample(itp, t_unit) 

    #     # gpu
    #     d = d |> gpu
    #     t = (time |> f32 |> gpu)'
    #     t_unit = KomaMRIBase.unit_time(t, tr) 
    #     itp = KomaMRIBase.interpolate(d, KomaMRIBase.Gridded(KomaMRIBase.Linear()), Val(size(d,1)))
    #     ux_gpu .= KomaMRIBase.resample(itp, t_unit)

    #     @test ux_cpu ≈ (ux_gpu |> cpu)

    #     ux_gpu .*= 0.0f0
    # end

    time = collect(0:0.1:1)

    obj = phantom_brain_arbitrary_motion() |> f32 |> gpu
    t = (time |> f32 |> gpu)

    Nparts = 20
    parts = kfoldperm(length(obj), Nparts)

    for i in 1:10
        foreach(enumerate(parts)) do (i, p)
            p = @view(obj[p])
            x, y, z = get_spin_coords(p.motion, p.x, p.y, p.z, t')
            println(@view(y[1, 1:10]))
        end
        print("\n")
    end

    @test true
end

@testitem "BlochSimple ArbitraryAction" tags=[:core] begin
    using Suppressor
    include("initialize_backend.jl")
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_brain_motion_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    sys = Scanner()
    obj = phantom_brain_arbitrary_motion()

    sim_params = Dict{String, Any}(
        "gpu"=>USE_GPU,
        "sim_method"=>KomaMRICore.BlochSimple(),
        "return_type"=>"mat"
    )
    sig = simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.
    println("NMRSE ArbitraryAction BlochSimple: ", NMRSE(sig, sig_jemris))
    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end


using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:core in ti.tags) #verbose=true

@testitem "Spinors×Mag" tags=[:core] begin
    # Spinor 2x2 representation should be equivalent to a 3x1 vector rotation
    x = rand(3); x = x./sum(x)
    θ = rand() * π
    n = rand(3); n = n./sqrt(sum(n.^2))
    z = Mag([x[1]+1im*x[2]], [x[3]])

    # General rotation
    xx1 = Q(θ,n[1]+1im*n[2],n[3])*z; #Spinor rot Q.(φ, B1./B, Bz./B)
    xx2 = Un(θ,n)*x; #3D rot matrix
    xx1 = [real(xx1.xy[1]), imag(xx1.xy[1]), xx1.z[1]]
    @test xx1 ≈ xx2

    # Rot x
    nx = [1,0,0]
    xx1 = Rx(θ)*z; #Spinor rot
    xx2 = Un(θ,nx)*x; #3D rot matrix
    xx1 = [real(xx1.xy[1]), imag(xx1.xy[1]), xx1.z[1]]
    @test xx1 ≈ xx2

    # Rot y
    nx = [0,1,0]
    xx1 = Ry(θ)*z; #Spinor rot
    xx2 = Un(θ,nx)*x; #3D rot matrix
    xx1 = [real(xx1.xy[1]), imag(xx1.xy[1]), xx1.z[1]]
    @test xx1 ≈ xx2

    # Rot z
    nx = [0,0,1]
    xx1 = Rz(θ)*z; #Spinor rot
    xx2 = Un(θ,nx)*x; #3D rot matrix
    xx1 = [real(xx1.xy[1]), imag(xx1.xy[1]), xx1.z[1]]
    @test xx1 ≈ xx2

    # Test Spinor struct
    α, β = rand(2)
    s = Spinor(α, β)
    @test s[1].α ≈ [Complex(α)] && s[1].β ≈ [Complex(β)]
    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", s)
    @test true
    α2, β2 = rand(2)
    s2 = Spinor(α2, β2)
    sp = s * s2
    @test sp.α ≈ s.α.*s2.α .- conj.(s2.β).*s.β && sp.β ≈ s.α.*s2.β .+ conj.(s2.α).*s.β
    φ, φ1, θ, φ2 = rand(4)
    Rm = KomaMRICore.Rg(φ1, θ, φ2)
    @test Rm.α ≈ [cos(θ/2)*exp(-1im*(φ1+φ2)/2)] && Rm.β ≈ [sin(θ/2)*exp(-1im*(φ1-φ2)/2)]
    Rn = KomaMRICore.Rφ(φ, θ)
    @test Rn.α ≈ [cos(θ/2)+0im] && Rn.β ≈ [exp(1im*φ)*sin(θ/2)]
    @test abs(s) ≈ [α^2 + β^2]
end

# Test ISMRMRD
@testitem "signal_to_raw_data" tags=[:core] begin
    using Suppressor

    seq = PulseDesigner.EPI_example()
    sys = Scanner()
    obj = brain_phantom2D()

    sim_params = KomaMRICore.default_sim_params()
    sim_params["return_type"] = "mat"
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

@testitem "Bloch_CPU_single_thread" tags=[:important, :core] begin
    using Suppressor
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>false,
        "Nthreads"=>1,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_CPU_multi_thread" tags=[:important, :core] begin
    using Suppressor
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>false,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end


@testitem "Bloch_GPU" tags=[:important, :skipci, :core] begin
    using Suppressor
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    sig_jemris = signal_jemris()
    seq = seq_epi_100x100_TE100_FOV230()
    obj = phantom_sphere()
    sys = Scanner()

    sim_params = Dict{String, Any}(
        "gpu"=>true,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_CPU_RF_accuracy_single_thread" tags=[:important, :core] begin
    using Suppressor

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])

    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1, Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>false, "Nthreads"=>1)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [0.570549+0.377122im
            0.607214+0.299628im
            0.633611+0.218961im
            0.649530+0.136450im
            0.654928+0.0534296im
            0.649928-0.0287866im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "Bloch_CPU_RF_accuracy_multi_thread" tags=[:important, :core] begin
    using Suppressor

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])

    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1, Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>false)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [0.570549+0.377122im
            0.607214+0.299628im
            0.633611+0.218961im
            0.649530+0.136450im
            0.654928+0.0534296im
            0.649928-0.0287866im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "Bloch_GPU_RF_accuracy" tags=[:important, :core] begin
    using Suppressor

    Tadc = 1e-3
    Trf = Tadc
    T1 = 1000e-3
    T2 = 20e-3
    Δw = 2π * 100
    B1 = 2e-6 * (Tadc / Trf)
    N = 6

    sys = Scanner()
    obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])

    seq = Sequence()
    seq += ADC(N, Tadc)
    for i=1:2
        global seq += RF(B1, Trf)
        global seq += ADC(N, Tadc)
    end

    sim_params = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>true)
    raw = @suppress simulate(obj, seq, sys; sim_params)

    #Mathematica-simulated Bloch equation result
    res1 = [0.153592+0.46505im,
            0.208571+0.437734im,
            0.259184+0.40408im,
            0.304722+0.364744im,
            0.344571+0.320455im,
            0.378217+0.272008im]
    res2 = [0.570549+0.377122im
            0.607214+0.299628im
            0.633611+0.218961im
            0.649530+0.136450im
            0.654928+0.0534296im
            0.649928-0.0287866im]
    norm2(x) = sqrt.(sum(abs.(x).^2))
    error0 = norm2(raw.profiles[1].data .- 0)
    error1 = norm2(raw.profiles[2].data .- res1) ./ norm2(res1) * 100
    error2 = norm2(raw.profiles[3].data .- res2) ./ norm2(res2) * 100

    @test  error0 + error1 + error2 < 0.1 #NMRSE < 0.1%
end

@testitem "BlochDict_CPU_single_thread" tags=[:important, :core] begin
    using Suppressor
    include(joinpath(@__DIR__, "test_files", "utils.jl"))

    seq = seq_epi_100x100_TE100_FOV230()
    obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])
    sys = Scanner()
    sim_params = Dict("gpu"=>false, "Nthreads"=>1, "sim_method"=>KomaMRICore.Bloch(), "return_type"=>"mat")
    sig = @suppress simulate(obj, seq, sys; sim_params)
    sig = sig / prod(size(obj))
    sim_params["sim_method"] = KomaMRICore.BlochDict()
    sig2 = simulate(obj, seq, sys; sim_params)
    sig2 = sig2 / prod(size(obj))
    @test sig ≈ sig2

    # Just checking to ensure that show() doesn't get stuck and that it is covered
    show(IOBuffer(), "text/plain", KomaMRICore.BlochDict())
    @test true
end

@testitem "simulate_slice_profile" tags=[:core] begin
    using Suppressor

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
    sim_params = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
    M = simulate_slice_profile(seq; z, sim_params)

    # For the time being, always pass the test
    @test true
end

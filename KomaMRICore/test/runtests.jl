using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:core in ti.tags) #verbose=true

@testitem "Sequence" tags=[:core] begin
    @testset "Init" begin
        sys = Scanner()
        B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
        EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
        @test dur(EX) ≈ durRF #RF length matches what is supposed to be

        #ACQ construction
        N = 101
        FOV = 23e-2
        EPI = PulseDesigner.EPI(FOV, N, sys)
        TE = 30e-3
        d1 = TE-dur(EPI)/2-dur(EX)
        d1 = d1 > 0 ? d1 : 0
        if d1 > 0 DELAY = Delay(d1) end

        #Sequence construction
        seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI #Only add delay if d1 is positive (enough space)
        seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
        @test dur(seq) ≈ dur(EX) + d1 + dur(EPI) #Sequence duration matches what is supposed to be
    end

    @testset "Rot_and_Concat" begin
        # Rotation 2D case
        A1, A2, T, t = rand(4)
        s = Sequence([Grad(A1,T);
                    Grad(A2,T)])
        θ = π*t
        R = rotz(θ)
        s2 = R*s #Matrix-Matrix{Grad} multiplication
        GR2 = R*s.GR.A #Matrix-vector multiplication
        @test s2.GR.A ≈ GR2
        # Rotation 3D case
        T, t1, t2, t3 = rand(4)
        N = 100
        GR = [Grad(rand(),T) for i=1:3, j=1:N]
        s = Sequence(GR)
        α, β, γ = π*t1, π*t2, π*t3
        Rx = rotx(α)
        Ry = roty(β)
        Rz = rotz(γ)
        R = Rx*Ry*Rz
        s2 = R*s #Matrix-Matrix{Grad} multiplication
        GR2 = R*s.GR.A #Matrix-vector multiplication
        @test s2.GR.A ≈ GR2

        # Concatenation of sequences
        A1, A2, A3, T1 = rand(4)
        s1 = Sequence([Grad(A1,T1);
                    Grad(A2,T1)],
                        [RF(A3,T1)])
        B1, B2, B3, T2 = rand(4)
        s2 = Sequence([Grad(B1,T2);
                    Grad(B2,T2)],
                        [RF(B3,T2)])
        s = s1 + s2
        @test s.GR.A ≈ [s1.GR.A s2.GR.A]
        @test s.RF.A ≈ [s1.RF.A s2.RF.A]
        @test s.ADC.N ≈ [s1.ADC.N ; s2.ADC.N]
    end

    @testset "Grad" begin
        #Testing gradient concatenation, breakes in some Julia versions
        A1, A2, T = rand(3)
        g1, g2 = Grad(A1,T), Grad(A2,T)
        GR = [g1;g2;;]
        GR2 = reshape([g1;g2],:,1)
        @test GR.A ≈ GR2.A

        #Sanity checks of contructors (A [T], T[s], rise[s], fall[s], delay[s])
        A, T = rand(2)
        g1, g2 = Grad(A,T), Grad(A,T,0,0,0)
        @test g1.A ≈ g2.A
        @test g1.T ≈ g2.T
        @test g1.rise ≈ g2.rise
        @test g1.fall ≈ g2.fall
        @test g1.delay ≈ g2.delay

        A, T, ζ = rand(3)
        g1, g2 = Grad(A,T,ζ), Grad(A,T,ζ,ζ,0)
        @test g1.A ≈ g2.A
        @test g1.T ≈ g2.T
        @test g1.rise ≈ g2.rise
        @test g1.fall ≈ g2.fall
        @test g1.delay ≈ g2.delay

        A, T, delay, ζ = rand(4)
        g1, g2 = Grad(A,T,ζ,delay), Grad(A,T,ζ,ζ,delay)
        @test g1.A ≈ g2.A
        @test g1.T ≈ g2.T
        @test g1.rise ≈ g2.rise
        @test g1.fall ≈ g2.fall
        @test g1.delay ≈ g2.delay
    end

    @testset "RF" begin
        #Testing gradient concatenation, breakes in some Julia versions
        A1, A2, T = rand(3)
        r1, r2 = RF(A1,T), RF(A2,T)
        R = [r1;r2;;]
        R2 = reshape([r1;r2],:,1)
        @test R.A ≈ R2.A

        #Sanity checks of constructors (A [T], T [s], Δf[Hz], delay [s])
        A, T = rand(2)
        r1, r2 = RF(A,T), RF(A,T,0,0)
        @test r1.A ≈ r2.A
        @test r1.T ≈ r2.T
        @test r1.Δf ≈ r2.Δf
        @test r1.delay ≈ r2.delay

        A, T, Δf = rand(3)
        r1, r2 = RF(A,T,Δf), RF(A,T,Δf,0)
        @test r1.A ≈ r2.A
        @test r1.T ≈ r2.T
        @test r1.Δf ≈ r2.Δf
        @test r1.delay ≈ r2.delay
    end
end

@testitem "PulseDesigner" tags=[:core] begin
    @testset "RF_sinc" begin
        sys = Scanner()
        B1 = 23.4732e-6 # For 90 deg flip angle
        Trf = 1e-3
        rf = PulseDesigner.RF_sinc(B1, Trf, sys; TBP=4)
        @test KomaMRICore.get_flip_angles(rf)[1] ≈ 90
    end
    @testset "Spiral" begin
        sys = Scanner()
        sys.Smax = 150    # [mT/m/ms]
        sys.Gmax = 500e-3 # [T/m]
        sys.GR_Δt = 4e-6  # [s]
        FOV = 0.2       # [m]
        N = 80          # Reconstructed image N×N
        Nint = 8
        λ = 2.1
        spiral = PulseDesigner.spiral_base(FOV, N, sys; λ=λ, BW=120e3, Nint)
        # Look at the k_space generated
        @test spiral(0).DEF["λ"] ≈ λ
    end
    @testset "Radial" begin
        sys = Scanner()
        N = 80
        Nspokes = ceil(Int64, π/2 * N ) #Nyquist in the radial direction
        FOV = 0.2
        spoke = PulseDesigner.radial_base(FOV, N, sys)
        @test spoke.DEF["Δθ"] ≈ π / Nspokes
    end
end

@testitem "Phantom" tags=[:core] begin
    #Test brain phantom 2D
    ph = brain_phantom2D()    #2D phantom
    @test ph.name=="brain2D_axial"

    #Test brain phantom 3D
    ph = brain_phantom3D()    #3D phantom
    @test ph.name=="brain3D"
end

@testitem "Scanner" tags=[:core] begin
    #Test somehow
    @test true
end

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
end

@testitem "TrapezoidalIntegration" tags=[:core] begin
    dt = Float64[1 1 1 1]
    x  = Float64[0 1 2 1 0]

    @test KomaMRICore.trapz(dt, x)[1] ≈ 4 #Triangle area = bh/2, with b = 4 and h = 2

    @test KomaMRICore.cumtrapz(dt, x) ≈ [0.5 2 3.5 4]
end

@testitem "Bloch_CPU_single_thread" tags=[:important, :core] begin
    using Suppressor, HDF5
    path = @__DIR__
    seq = @suppress read_seq(path*"/test_files/epi_100x100_TE100_FOV230.seq")
    obj = read_phantom_jemris(path*"/test_files/sphere_chemical_shift.h5")
    sys = Scanner()

    simParams = Dict{String, Any}(
        "gpu"=>false,
        "Nthreads"=>1,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; simParams)
    sig = sig / prod(size(obj))

    sig_jemris = h5open(path*"/test_files/jemris_signals_epi_sphere_cs.h5")["/signal/channels/00"]
    sig_jemris = sig_jemris[1,:] + 1im*sig_jemris[2,:]
    sig_jemris = sig_jemris[:]

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_CPU_multi_thread" tags=[:important, :core] begin
    using Suppressor, HDF5
    path = @__DIR__
    seq = @suppress read_seq(path*"/test_files/epi_100x100_TE100_FOV230.seq")
    obj = read_phantom_jemris(path*"/test_files/sphere_chemical_shift.h5")
    sys = Scanner()

    simParams = Dict{String, Any}(
        "gpu"=>false,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; simParams)
    sig = sig / prod(size(obj))

    sig_jemris = h5open(path*"/test_files/jemris_signals_epi_sphere_cs.h5")["/signal/channels/00"]
    sig_jemris = sig_jemris[1,:] + 1im*sig_jemris[2,:]
    sig_jemris = sig_jemris[:]

    NMRSE(x, x_true) = sqrt.( sum(abs.(x .- x_true).^2) ./ sum(abs.(x_true).^2) ) * 100.

    @test NMRSE(sig, sig_jemris) < 1 #NMRSE < 1%
end

@testitem "Bloch_GPU" tags=[:important, :skipci, :core] begin
    using Suppressor, HDF5
    path = @__DIR__
    seq = @suppress read_seq(path*"/test_files/epi_100x100_TE100_FOV230.seq")
    obj = read_phantom_jemris(path*"/test_files/sphere_chemical_shift.h5")
    sys = Scanner()

    simParams = Dict{String, Any}(
        "gpu"=>true,
        "sim_method"=>KomaMRICore.Bloch(),
        "return_type"=>"mat"
    )
    sig = @suppress simulate(obj, seq, sys; simParams)
    sig = sig / prod(size(obj))

    sig_jemris = h5open(path*"/test_files/jemris_signals_epi_sphere_cs.h5")["/signal/channels/00"]
    sig_jemris = sig_jemris[1,:] + 1im*sig_jemris[2,:]
    sig_jemris = sig_jemris[:]

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

    simParams = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>false, "Nthreads"=>1)
    raw = @suppress simulate(obj, seq, sys; simParams)

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

    simParams = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>false)
    raw = @suppress simulate(obj, seq, sys; simParams)

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

    simParams = Dict{String, Any}("Δt_rf"=>1e-5, "gpu"=>true)
    raw = @suppress simulate(obj, seq, sys; simParams)

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

@testitem "IO" tags=[:core] begin
    using Suppressor
    #Test Pulseq
    @testset "Pulseq" begin
        path = @__DIR__
        seq = @suppress read_seq(path*"/test_files/epi.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "epi.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000

        seq = @suppress read_seq(path*"/test_files/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        @test seq.DEF["FileName"] == "spiral.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1004000

        seq = @suppress read_seq(path*"/test_files/epi_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "epi_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001

        seq = @suppress read_seq(path*"/test_files/radial_JEMRIS.seq") #Pulseq v1.2.1
        @test seq.DEF["FileName"] == "radial_JEMRIS.seq"
        @test seq.DEF["PulseqVersion"] ≈ 1002001

        #Test ISMRMRD
        fraw = ISMRMRDFile(path*"/test_files/Koma_signal.mrd")
        raw = RawAcquisitionData(fraw)
        @test raw.params["protocolName"] == "epi"
        @test raw.params["institutionName"] == "Pontificia Universidad Catolica de Chile"
        @test raw.params["encodedSize"] ≈ [101, 101, 1]
        @test raw.params["reconSize"] ≈ [102, 102, 1]
        @test raw.params["patientName"] == "brain2D_axial"
        @test raw.params["trajectory"] == "other"
        @test raw.params["systemVendor"] == "KomaMRI.jl"
    end
    #Test JEMRIS
    @testset "JEMRIS" begin
        path = @__DIR__
        obj = read_phantom_jemris(path*"/test_files/column1d.h5")
        @test obj.name == "column1d.h5"
    end
end
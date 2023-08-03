using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:core in ti.tags) #verbose=true

@testitem "Sequence" tags=[:core] begin
    using Suppressor
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
        A, T = 0.1, 1e-3
        grad = Grad(A, T)

        A, T = rand(2)
        g1, g2 = Grad(A,T), Grad(A,T,0,0,0)
        @test g1 ≈ g2

        A, T, ζ = rand(3)
        g1, g2 = Grad(A,T,ζ), Grad(A,T,ζ,ζ,0)
        @test g1 ≈ g2

        A, T, delay, ζ = rand(4)
        g1, g2 = Grad(A,T,ζ,delay), Grad(A,T,ζ,ζ,delay)
        @test g1 ≈ g2

        # Test construction with shape function
        T, N = 1e-3, 100
        f = t -> sin(π*t / T)
        gradw = Grad(f, T, N)
        @test gradw.A ≈ f.(range(0, T; length=N))

        # Test Grad operations
        α = 3
        gradt = α * grad
        @test size(grad, 1) == 1
        @test gradt.A ≈ α * grad.A
        gradt = grad * α
        @test gradt.A ≈ α * grad.A
        gradt = grad / α
        @test gradt.A ≈ grad.A / α
        grads = grad + gradt
        @test grads.A ≈ grad.A + gradt.A
        A1, A2, A3 = 0.1, 0.2, 0.3
        v1 = [Grad(A1,T); Grad(A2,T); Grad(A3,T)]
        v2 = [Grad(A2,T); Grad(A3,T); Grad(A1,T)]
        v3 = v1 + v2
        @test [v3[i].A for i=1:length(v3)] ≈ [v1[i].A + v2[i].A for i=1:length(v1)]
        gradr = grad - gradt
        @test gradr.A ≈ grad.A - gradt.A
        gradt = -grad
        @test gradt.A ≈ -grad.A
        vc = vcat(v1, v2)
        @test [vc[1,j].A for j=1:length(v1)] ≈ [v1[i].A for i=1:length(v1)]
        @test [vc[2,j].A for j=1:length(v2)] ≈ [v2[i].A for i=1:length(v2)]
        vc = vcat(v1, v2, v3)
        @test [vc[1,j].A for j=1:length(v1)] ≈ [v1[i].A for i=1:length(v1)]
        @test [vc[2,j].A for j=1:length(v2)] ≈ [v2[i].A for i=1:length(v2)]
        @test [vc[3,j].A for j=1:length(v3)] ≈ [v3[i].A for i=1:length(v3)]
        delay, rise, T, fall = 1e-6, 2e-6, 10e-3, 3e-6
        gr = Grad(A, T, rise, fall, delay)
        @test dur(gr) ≈ delay + rise + T + fall
        T1, T2, T3 = 1e-3, 2e-3, 3e-3
        vt = [Grad(A1,T1); Grad(A2,T2); Grad(A3,T3)]
        @test dur(vt) ≈ [maximum([T1, T2, T3])]

        # Test Grad output message
        io = IOBuffer()
        show(io, "text/plain", grad)
        @test occursin("Grad(", String(take!(io)))
        show(io, "text/plain", gr)
        @test occursin("Grad(", String(take!(io)))

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
        @test r1 ≈ r2

        A, T, Δf = rand(3)
        r1, r2 = RF(A,T,Δf), RF(A,T,Δf,0)
        @test r1 ≈ r2

        # Test RF output message
        io = IOBuffer()
        show(io, "text/plain", r1)
        @test occursin("RF(", String(take!(io)))

        # Test Grad operations
        B1x, B1y, T = rand(3)
        A = B1x + im*B1y
        α = Complex(rand())
        rf = RF(A, T)
        rft = α * rf
        @test size(rf, 1) == 1
        @test rft.A ≈ α * rf.A
        @test dur(rf) ≈ rf.T
        B1x, B1y, B2x, B2y, B3x, B3y, T1, T2, T3 = rand(9)
        rf1, rf2, rf3 = RF(B1x + im*B1y, T1), RF(B1x + im*B1y, T2), RF(B3x + im*B3y, T3)
        rv = [rf1; rf2; rf3]
        @test dur(rv) ≈ sum(dur.(rv))

    end

    @testset "Delay" begin

        # Test delay construction
        T = 1e-3
        delay = Delay(T)
        @test delay.T ≈ T

        # Test delay construction error for negative values
        err = Nothing
        try Delay(-T) catch err end
        @test err isa ErrorException

        # Test delay output message
        io = IOBuffer()
        show(io, "text/plain", delay)
        @test occursin("Delay(", String(take!(io)))

        # Test addition of a delay to a sequence
        seq = Sequence()
        ds = delay + seq
        @test dur(ds[1]) ≈ delay.T && dur(ds[2]) ≈ .0
        sd = seq + delay
        @test dur(sd[1]) ≈ .0 && dur(sd[2]) ≈ delay.T

    end
    @testset "ADC" begin

        # Test ADC construction
        N, T, delay, Δf, ϕ  = 64, 1e-3, 2e-3, 1e-6, .25*π
        adc = ADC(N, T, delay, Δf, ϕ)

        adc1, adc2 = ADC(N, T), ADC(N,T,0,0,0)
        @test adc1 ≈ adc2

        adc1, adc2 = ADC(N, T, delay), ADC(N, T, delay, 0, 0)
        @test adc1 ≈ adc2

        adc1, adc2 = ADC(N, T, delay, Δf, ϕ), ADC(N, T, delay, Δf, ϕ)
        @test adc1 ≈ adc2

        # Test ADC construction errors for negative values
        err = Nothing
        try ADC(N, -T) catch err end
        @test err isa ErrorException
        try ADC(N, -T,  delay) catch err end
        @test err isa ErrorException
        try ADC(N,  T, -delay) catch err end
        @test err isa ErrorException
        try ADC(N, -T, -delay) catch err end
        @test err isa ErrorException
        try ADC(N, -T,  delay, Δf, ϕ) catch err end
        @test err isa ErrorException
        try ADC(N,  T, -delay, Δf, ϕ) catch err end
        @test err isa ErrorException
        try ADC(N, -T, -delay, Δf, ϕ) catch err end
        @test err isa ErrorException

        # Test ADC getproperties
        Nb, Tb, delayb, Δfb, ϕb  = 128, 2e-3, 4e-3, 2e-6, .125*π
        adb = ADC(Nb, Tb, delayb, Δfb, ϕb)
        adcs = [adc, adb]
        @test adcs.N ≈ [adc.N, adb.N] && adcs.T ≈ [adc.T, adb.T] && adcs.delay ≈ [adc.delay, adb.delay]
        @test adcs.Δf ≈ [adc.Δf, adb.Δf] && adcs.ϕ ≈ [adc.ϕ, adb.ϕ] && adcs.dur ≈ [adc.T + adc.delay, adb.T + adb.delay]

    end

    @testset "DiscreteSequence" begin
        path = @__DIR__
        seq = @suppress read_seq(path*"/test_files/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        simParams = KomaMRICore.default_sim_params()
        t, Δt = KomaMRICore.get_uniform_times(seq, simParams["Δt"]; Δt_rf=simParams["Δt_rf"])
        seqd = KomaMRICore.discretize(seq)
        i1, i2 = rand(1:Int(floor(0.5*length(seqd)))), rand(Int(ceil(0.5*length(seqd))):length(seqd))
        @test seqd[i1].t ≈ [t[i1]]
        @test seqd[i1:i2-1].t ≈ t[i1:i2]

        T, N = 1, 4
        seq = RF(1, 1)
        seq += Sequence([Grad(1, 1)])
        seq += ADC(N, 1)
        simParams = KomaMRICore.default_sim_params()
        simParams["Δt"], simParams["Δt_rf"] = T/N, T/N
        seqd = KomaMRICore.discretize(seq; simParams)
        i = Int(floor(length(seqd) / 3))
        @test is_RF_on(seq[1]) == is_RF_on(seqd[1*i:1*i+1]) && is_GR_on(seq[1]) == is_GR_on(seqd[1*i:1*i+1]) && is_ADC_on(seq[1]) == is_ADC_on(seqd[1*i:1*i+1])
        @test is_RF_on(seq[2]) == is_RF_on(seqd[2*i:2*i+1]) && is_GR_on(seq[2]) == is_GR_on(seqd[2*i:2*i+1]) && is_ADC_on(seq[2]) == is_ADC_on(seqd[2*i:2*i+1])
        @test is_RF_on(seq[3]) == is_RF_on(seqd[3*i:3*i+1]) && is_GR_on(seq[3]) == is_GR_on(seqd[3*i:3*i+1]) && is_ADC_on(seq[3]) == is_ADC_on(seqd[3*i:3*i+1])
        @test KomaMRICore.is_GR_off(seqd) ==  !KomaMRICore.is_GR_on(seqd)
        @test KomaMRICore.is_RF_off(seqd) ==  !KomaMRICore.is_RF_on(seqd)
        @test KomaMRICore.is_ADC_off(seqd) == !KomaMRICore.is_ADC_on(seqd)
    end

    @testset "SequenceFunctions" begin
        path = @__DIR__
        seq = @suppress read_seq(path*"/test_files/spiral.seq") #Pulseq v1.4.0, RF arbitrary
        t, Δt = KomaMRICore.get_uniform_times(seq, 1)
        t_adc =  KomaMRICore.get_adc_sampling_times(seq)
        M2, M2_adc = KomaMRICore.get_slew_rate(seq)
        Gx, Gy, Gz = KomaMRICore.get_grads(seq, t)
        Gmx, Gmy, Gmz = KomaMRICore.get_grads(seq, reshape(t, 1, :))
        @test reshape(Gmx, :, 1) ≈ Gx && reshape(Gmy, :, 1) ≈ Gy && reshape(Gmz, :, 1) ≈ Gz
        @test is_ADC_on(seq) == is_ADC_on(seq, t)
        @test is_RF_on(seq) == is_RF_on(seq, t)
        @test KomaMRICore.is_Delay(seq) == !(is_GR_on(seq) || is_RF_on(seq) || is_ADC_on(seq))
        @test size(M2, 1) == length(Δt) && size(M2_adc, 1) == length(t_adc)
        io = IOBuffer()
        show(io, "text/plain", seq)
        @test occursin("Sequence[", String(take!(io)))

        α = rand()
        c = α + im*rand()
        x = seq
        y = @suppress read_seq(path*"/test_files/epi.seq") #Pulseq v1.4.0, RF arbitrary
        z = x + y
        @test z.GR.A ≈ [x.GR  y.GR].A && z.RF.A ≈ [x.RF y.RF].A && z.ADC.N ≈ [x.ADC; y.ADC].N
        z = x - y
        @test z.GR.A ≈ [x.GR  -y.GR].A
        z = -x
        @test z.GR.A ≈ -x.GR.A
        z = x * α
        @test z.GR.A ≈ α*x.GR.A
        z = α * x
        @test z.GR.A ≈ α*x.GR.A
        z = x * c
        @test z.RF.A ≈ c*x.RF.A
        z = c * x
        @test z.RF.A ≈ c*x.RF.A
        z = x / α
        @test z.GR.A ≈ x.GR.A/α
        @test size(y) == size(y.GR[1,:])
        z = x + x.GR[3,1]
        @test z.GR.A[1, end] ≈ x.GR[3,1].A
        z = x.GR[3,1] + x
        @test z.GR.A[1, 1] ≈ x.GR[3,1].A
        z = x + x.RF[1,1]
        @test z.RF.A[1, end] ≈ x.RF[1,1].A
        z = x.RF[1,1] + x
        @test z.RF.A[1, 1] ≈ x.RF[1,1].A
        z = x + x.ADC[3,1]
        @test z.ADC.N[end] ≈ x.ADC[3,1].N
        z = x.ADC[3,1] + x
        @test z.ADC.N[1] ≈ x.ADC[3,1].N
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

    # Test phantom struct creation
    name = "Bulks"
    x = [-2e-3; -1e-3; 0.; 1e-3; 2e-3]
    y = [-4e-3; -2e-3; 0.; 2e-3; 4e-3]
    z = [-6e-3; -3e-3; 0.; 3e-3; 6e-3]
    ρ = [.2; .4; .6; .8; 1.]
    T1 = [.9; .9; .5; .25; .4]
    T2 = [.09; .05; .04; .07; .005]
    T2s = [.1; .06; .05; .08; .015]
    Δw = [-2e-6; -1e-6; 0.; 1e-6; 2e-6]
    Dλ1 = [-4e-6; -2e-6; 0.; 2e-6; 4e-6]
    Dλ2 = [-6e-6; -3e-6; 0.; 3e-6; 6e-6]
    Dθ = [-8e-6; -4e-6; 0.; 4e-6; 8e-6]
    u = (x,y,z,t)->0
    obj = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, ux=u, uy=u, uz=u)
    obj2 = Phantom(name=name, x=x, y=y, z=z, ρ=ρ, T1=T1, T2=T2, T2s=T2s, Δw=Δw, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, ux=u, uy=u, uz=u)
    @test obj ≈ obj2

    # Test size and length definitions of a phantom
    @test size(obj) == size(ρ)
    #@test length(obj) == length(ρ)

    # Test phantom comparison
    ue = (x,y,z,t)->1
    obe = Phantom(name, x, y, z, ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, ue, ue, ue)
    @test obj ≈ obe

    # Test phantom subset
    rng = 1:2:5
    ug = (x,y,z,t)->-1
    obg = Phantom(name, x[rng], y[rng], z[rng], ρ[rng], T1[rng], T2[rng], T2s[rng], Δw[rng],
                    Dλ1[rng], Dλ2[rng], Dθ[rng], ug, ug, ug)
    @test obj[rng] ≈ obg
    @test @view(obj[rng]) ≈ obg

    # Test addition of phantoms
    ua = (x,y,z,t)->2
    oba = Phantom(name, [x; x[rng]], [y; y[rng]], [z; z[rng]], [ρ; ρ[rng]],
                    [T1; T1[rng]], [T2; T2[rng]], [T2s; T2s[rng]], [Δw; Δw[rng]],
                    [Dλ1; Dλ1[rng]], [Dλ2; Dλ2[rng]], [Dθ; Dθ[rng]], ua, ua, ua)
    @test obj + obg ≈ oba

    # Test scalar multiplication of a phantom
    c = 7.
    obc = Phantom(name, x, y, z, c*ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, u, u, u)
    @test c*obj ≈ obc

    #Test brain phantom 2D
    ph = brain_phantom2D()
    @test ph.name=="brain2D_axial"

    #Test brain phantom 3D
    ph = brain_phantom3D()
    @test ph.name=="brain3D"
end

@testitem "Scanner" tags=[:core] begin
    B0, B1, Gmax, Smax = 1.5, 10e-6, 60e-3, 500
    ADC_Δt, seq_Δt, GR_Δt, RF_Δt = 2e-6, 1e-5, 1e-5, 1e-6
    RF_ring_down_T, RF_dead_time_T, ADC_dead_time_T = 20e-6, 100e-6, 10e-6
    sys = Scanner(B0, B1, Gmax, Smax, ADC_Δt, seq_Δt, GR_Δt, RF_Δt, RF_ring_down_T, RF_dead_time_T, ADC_dead_time_T)
    @test sys.B0 ≈ B0 && sys.B1 ≈ B1 && sys.Gmax ≈ Gmax && sys.Smax ≈ Smax
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

    # Test Spinor struct
    α, β = rand(2)
    s = Spinor(α, β)
    @test s[1].α ≈ [Complex(α)] && s[1].β ≈ [Complex(β)]
    io = IOBuffer()
    show(io, "text/plain", s)
    @test occursin("Spinor(", String(take!(io)))
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

@testitem "TimeStepCalculation" tags=[:core] begin
    ampRF = 1e-6
    durRF = 1e-3
    index_offset, number_rf_points = 2, 3
    rf_key_points = []
    rf = RF(ampRF, durRF)
    seq = Sequence()
    seq += rf
    append!(rf_key_points, [index_offset; index_offset + number_rf_points + 1])
    seq += Delay(durRF)
    seq += rf
    append!(rf_key_points, [rf_key_points[end] + 1; rf_key_points[end] + 1 + number_rf_points + 1])
    seq += Delay(durRF)
    seq += rf
    append!(rf_key_points, [rf_key_points[end] + 1; rf_key_points[end] + 1 + number_rf_points + 1])
    t, Δt = KomaMRICore.get_variable_times(seq; dt_rf=durRF)
    @test KomaMRICore.get_breaks_in_RF_key_points(seq, t) == rf_key_points
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

@testitem "BlochDict_CPU_single_thread" tags=[:important, :core] begin
    using Suppressor
    path = joinpath(@__DIR__, "test_files")
    seq = @suppress read_seq(joinpath(path, "epi_100x100_TE100_FOV230.seq"))
    obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])
    sys = Scanner()
    simParams = Dict("gpu"=>false, "Nthreads"=>1, "sim_method"=>KomaMRICore.Bloch(), "return_type"=>"mat")
    sig = @suppress simulate(obj, seq, sys; simParams)
    sig = sig / prod(size(obj))
    simParams["sim_method"] = KomaMRICore.BlochDict()
    sig2 = simulate(obj, seq, sys; simParams)
    sig2 = sig2 / prod(size(obj))
    @test sig ≈ sig2

    io = IOBuffer()
    show(io, "text/plain", KomaMRICore.BlochDict())
    @test occursin("BlochDict(", String(take!(io)))
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

        # Test signal_to_raw_data
        signal1 = Vector()
        for i=1:length(raw.profiles)
            signal1 = [signal1; raw.profiles[i].data]
        end
        rawmrd = signal_to_raw_data(signal1, seq)
        @test rawmrd.params["institutionName"] == raw.params["institutionName"]
        io = IOBuffer()
        show(io, "text/plain", rawmrd)
        @test occursin("RawAcquisitionData[", String(take!(io)))

    end
    #Test JEMRIS
    @testset "JEMRIS" begin
        path = @__DIR__
        obj = read_phantom_jemris(path*"/test_files/column1d.h5")
        @test obj.name == "column1d.h5"
    end
    #Test JEMRIS
    @testset "MRiLab" begin
        path = @__DIR__
        filename = path * "/test_files/brain_mrilab.mat"
        FRange_filename = path * "/test_files/FRange.mat" #Slab within slice thickness
        obj = read_phantom_MRiLab(filename; FRange_filename)
        @test obj.name == "brain_mrilab.mat"
    end
end

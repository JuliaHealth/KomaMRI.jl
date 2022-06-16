using KomaMRI
using Test

@testset "Sequence" begin
    @testset "Sequence_op" begin
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

@testset "Phantom" begin
    #Test somehow
end

@testset "Scanner" begin
    #Test somehow
end

@testset "Simulation" begin
    @testset "Spinors×Mag" begin
        # Spinor 2x2 representation should be equivalent to a 3x1 vector rotation
        x = rand(3); x = x./sum(x)
        θ = rand() * π
        n = rand(3); n = n./sqrt(sum(n.^2))
        z = Mag(x[1]+1im*x[2], x[3])
    
        # General rotation
        xx1 = Q(θ,n[1]+1im*n[2],n[3])*z; #Spinor rot Q.(φ, B1./B, Bz./B)
        xx2 = Un(θ,n)*x; #3D rot matrix
        xx1 = [real(xx1.xy), imag(xx1.xy), xx1.z]
        @test xx1 ≈ xx2
    
        # Rot x
        nx = [1,0,0]
        xx1 = Rx(θ)*z; #Spinor rot
        xx2 = Un(θ,nx)*x; #3D rot matrix
        xx1 = [real(xx1.xy), imag(xx1.xy), xx1.z]
        @test xx1 ≈ xx2
    
        # Rot y
        nx = [0,1,0]
        xx1 = Ry(θ)*z; #Spinor rot
        xx2 = Un(θ,nx)*x; #3D rot matrix
        xx1 = [real(xx1.xy), imag(xx1.xy), xx1.z]
        @test xx1 ≈ xx2
    
        # Rot z
        nx = [0,0,1]
        xx1 = Rz(θ)*z; #Spinor rot
        xx2 = Un(θ,nx)*x; #3D rot matrix
        xx1 = [real(xx1.xy), imag(xx1.xy), xx1.z]
        @test xx1 ≈ xx2
    end
    #Test somehow
end

@testset "Recon" begin
    #Test
end

@testset "IO" begin
    #Test Pulseq

    #Test ISMRMRD

    #Test JEMRIS
end

#GUI tests
@testset "GUI" begin
    @testset "GUI_phantom" begin
        ph = brain_phantom2D()   #2D phantom
        plot_phantom_map(ph, :ρ) #Plotting the phantom's rho map
        @test true               #If the previous line fails the test will fail

        plot_phantom_map(ph, :T1) #Plotting the phantom's rho map
        @test true               #If the previous line fails the test will fail

        plot_phantom_map(ph, :T2) #Plotting the phantom's rho map
        @test true               #If the previous line fails the test will fail
    end

    @testset "GUI_seq" begin
        #RF construction
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

        #Plot sequence
        plot_seq(seq)  #Plotting the sequence
        @test true          #If the previous line fails the test will fail

        #Plot k-space
        plot_kspace(seq)    #Plotting the k-space
        @test true          #If the previous line fails the test will fail

        #Plot M0
        plot_M0(seq)        #Plotting the M0
        @test true          #If the previous line fails the test will fail
    end
end

nothing
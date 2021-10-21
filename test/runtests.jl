using MRIsim
using Test

@testset "Spinors×Magnetization" begin
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

@testset "Sequence" begin
    # Rotation 2D case
    A1, A2, T, t = rand(4)
    s = Sequence([Grad(A1,T);
                  Grad(A2,T)])
    θ = π*t
    R = MRIsim.rotz(θ)
    s2 = R*s
    @test s2.GR.A ≈ [A1*cos(θ)-A2*sin(θ); A1*sin(θ)+A2*cos(θ)]
    # Rotation 3D case
    T, t1, t2, t3 = rand(4)
    N = 100
    GR = [Grad(rand(),T) for i=1:3, j=1:N]
    s = Sequence(GR)
    α, β, γ = π*t1, π*t2, π*t3
    Rx = MRIsim.rotx(α)
    Ry = MRIsim.roty(β)
    Rz = MRIsim.rotz(γ)
    R = Rx*Ry*Rz
    s2 = R*s
    GR2 = R*GR.A
    @test s2.GR.A ≈ GR2
    # Concatenation
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
    @test s.DAC.N ≈ [s1.DAC.N ; s2.DAC.N]
end
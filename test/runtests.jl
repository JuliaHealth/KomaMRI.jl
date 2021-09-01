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

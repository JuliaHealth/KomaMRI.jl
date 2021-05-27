using MRIsim
using Test

@testset "Spinors" begin
    # Spinor 2x2 representation should be equivalent to a 3x1 vector rotation
    x = rand(3); x = x./sum(x)
    θ = rand() * π
    n = rand(3); n = n./sqrt(sum(n.^2))
    z = [x[1]+1im*x[2], x[3]]
    # General rotation
    xx1 = MRIsim.Q(θ,n)*z; #Spinor rot
    xx2 = MRIsim.Un(θ,n)*x; #3D rot matrix
    xx1 = [real(xx1[1]), imag(xx1[1]), xx1[2]]
    @test xx1 ≈ xx2
    # Rot x
    nx = [1,0,0]
    xx1 = MRIsim.Rx(θ)*z; #Spinor rot
    xx2 = MRIsim.Un(θ,nx)*x; #3D rot matrix
    xx1 = [real(xx1[1]), imag(xx1[1]), xx1[2]]
    @test xx1 ≈ xx2
    # Rot y
    nx = [0,1,0]
    xx1 = MRIsim.Ry(θ)*z; #Spinor rot
    xx2 = MRIsim.Un(θ,nx)*x; #3D rot matrix
    xx1 = [real(xx1[1]), imag(xx1[1]), xx1[2]]
    @test xx1 ≈ xx2
    # Rot z
    nx = [0,0,1]
    xx1 = MRIsim.Rz(θ)*z; #Spinor rot
    xx2 = MRIsim.Un(θ,nx)*x; #3D rot matrix
    xx1 = [real(xx1[1]), imag(xx1[1]), xx1[2]]
    @test xx1 ≈ xx2
end

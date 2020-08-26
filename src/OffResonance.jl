# B = μ (I-σI+1/3Χ) H
#FT^-1{1/}

using Plots
using LinearAlgebra
pyplot()

Dipole(x,y) = begin
    r = sqrt.(x.^2 .+ y.^2)
    θ = atan.(x,y) #angle with respect to the main magnetic field
    r < 1e-16 ? 0 : 1/(4*π)*(1 .+ 3*cos.(2*θ))./r.^3
end

B0 = 1 #From scanner object

"""Schenck JF. The role of magnetic susceptibility in magnetic resonance
imaging: MRI magnetic compatibility of the first and second kinds.
Med Phys. 1996;23(6):815-850. doi:10.1118/1.597854"""
Sphere(x,R,Δχ) = begin
    if norm(x) <= R #Inside sphere
        2Δχ/3*B0 #2/3
    else #Outside
        Δχ/3*B0*R^3*(2x[3]^2-x[1]^2-x[2]^2)/norm(x)^5
    end
end

"""Schenck JF. The role of magnetic susceptibility in magnetic resonance
imaging: MRI magnetic compatibility of the first and second kinds.
Med Phys. 1996;23(6):815-850. doi:10.1118/1.597854"""
CylinderZ(x,R,Δχ) = begin
    if norm(x) <= R #Inside sphere
        Δχ*B0 #1
    else #Outside
        0
    end
end

"""Schenck JF. The role of magnetic susceptibility in magnetic resonance
imaging: MRI magnetic compatibility of the first and second kinds.
Med Phys. 1996;23(6):81 5-850. doi:10.1118/1.597854"""
CylinderY(x,R,Δχ) = begin
    if norm(x) <= R #Inside sphere
        Δχ/2*B0 #1/2
    else #Outside
        Δχ/2*B0*R^2*(x[3]^2-x[1]^2)/(x[3]^2+x[1]^2)^2
    end
end

contourf([CylinderY([x,0,z],.1,.1) for z=-1:.001:1, x=-1:.001:1],aspect_ratio=:equal)
contourf([Sphere([x,0,z],.1,.1) for z=-1:.001:1, x=-1:.001:1],aspect_ratio=:equal)

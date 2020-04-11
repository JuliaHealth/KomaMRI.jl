using Plots
import Base.*
"""
	Spinor(α,β)

Based on "Introduction to the Shinnar-Le Roux algorithm", Patrick Le Roux (1995).

These rotations are represented by 2 X 2 unitary matrix (α,β∈C):

R = [α -β⃰;

-------- β   α⃰],

and |α|²+|β|² = 1. Equivalent to Spinor (α,β) with its Cayley-Klein parameters.
"""
struct Spinor
	α::Complex
	β::Complex
end
"""
Spinor multiplication identity.
	(α1,β1)×(α2,β2) = (α1 α2 - β2⃰ β1 , β2 α1 + α2⃰ β1)
"""
*(s1::Spinor,s2::Spinor) = begin
	Spinor(s1.α*s2.α - conj(s2.β)*s1.β,
		   s1.α*s2.β + conj(s2.α)*s1.β)
end
"""
Spinor × Magnetization

V = (x,y,z)

<=>

V = [z -X⃰;

-------- X  -z],	X = x + i y.

Then

(α,β)×(X,z) = (X+,z+),

with

X+ = 2α⃰βz+α⃰²X-β²X⃰

and

z+ = (|α|²-|β|²)z-α⃰β⃰X-αβX⃰ .
"""
*(s::Spinor,M::Array) = begin
	[2*conj(s.α)*s.β*M[2]+conj(s.α)^2*M[1]-s.β^2*conj(M[1]),
						   (abs(s.α)^2-abs(s.β)^2)*M[2]-conj(s.α)*conj(s.β)*M[1]-s.α*s.β*conj(M[1])]
end
#Rotation matrices
Rz(φ) = Spinor(exp(-im*φ/2), 0)
Ry(θ) = Spinor(cos(θ/2), sin(θ/2))
Rx(θ) = Spinor(cos(θ/2), -im*sin(θ/2)) #Rz(π/2)*Ry(θ)*Rz(-π/2)
Rg(φ1,θ,φ2) = Spinor(cos(θ/2)*exp(-im*(φ1+φ2)/2), sin(θ/2)*exp(-im*(φ1-φ2)/2))#Rz(φ2)*Ry(θ)*Rz(φ1)
Rφ(φ,θ) = Spinor(cos(θ/2),exp(-im*φ)*sin(θ/2))

"""
Pauly, J., Le Roux, P., Nishimura, D., & Macovski, A. (1991).
Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm (NMR imaging).
IEEE Transactions on Medical Imaging, 10(1), 53–65. doi:10.1109/42.75611

Spin-domain rotation matrix.

Rotation of φ, with respect to the axis of rotation n=(nx,ny,nz).

φ = -γ Δt √(|B1j|²+(G⋅x)²)

n =  γ Δt/|φ| (B1x, B1y, G⋅x)
"""
Q(φ,n::Array{Float64}) = begin
	Spinor(cos(φ/2)-im*n[3]*sin(φ/2),
			-im*(n[1]+im*n[2])*sin(φ/2))
end

# ## Example
# G = 30e-3
# γ = 2*π*42.5e6
# T = 3e-3; N = 200
# t = 0:T/(N-1):T; Δt = t[2]-t[1]
#
# B1 = 15e-6; B1e = B1*sinc.(8(t.-T/2)/T)
# φ = -γ*Δt*sqrt.(abs.(B1e).^2) ; φ[φ.==0] .= 1e-17
# n =  γ*Δt*[[B1e[i]/abs(φ[i]) 0 0] for i = 1:N]
# Qs = [Q(φ[i],n[i]) for i=1:N]
# Qt = *(Qs...) #Total rotation matrix
#
# Mxy, Mz = 0. + 0. *im, 1. #[0,0,1]
# M = [Mxy, Mz]
# Mt = Qt*M
# Mz = [(*(Qs[1:i]...)*M)[2] for i=2:N]
# Mxy = [(*(Qs[1:i]...)*M)[1] for i=2:N]
# #Sphere
# u = range(0,stop=2*π,length=100);
# v = range(0,stop=π,length=100);
# x = .99 * cos.(u) * sin.(v)';
# y = .99 * sin.(u) * sin.(v)';
# z = .99 * ones(100 ) * cos.(v)';
#
# plotlyjs()
# l = @layout [a ;  b]
# p1 = plot(B1e/B1,label="B1e")
# plot!(real(Mz),label="Mz")
# plot!(abs.(Mxy),label="Mxy")
# p2 = plot(real(Mxy),imag(Mxy),real(Mz),xlim=(-1,1),ylim=(-1,1),zlim=(-1,1),
# 	linewidth=5,colorbar=false,legend=false)
# surface!(x,y,z,color="white",alpha=.8)
# plot(p1,p2,layout=l)

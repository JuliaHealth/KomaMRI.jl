"""
	Spinor(α,β)

Based on "Introduction to the Shinnar-Le Roux algorithm", Patrick Le Roux (1995).

These rotations are represented by 2 X 2 unitary matrix (α,β ∈ C):

R = [α -β⃰;

-------- β   α⃰],

and |α|²+|β|² = 1. Equivalent to Spinor (α,β) with its Cayley-Klein parameters.
"""
struct Spinor
	α::Complex
	β::Complex
end
Base.show(io::IO,s::Spinor) = print(io, "(α = ",s.α,", β = ",s.β,")")

"""
Spinor multiplication identity.
	(α1,β1)×(α2,β2) = (α1 α2 - β2⃰ β1 , β2 α1 + α2⃰ β1)
"""
*(s1::Spinor,s2::Spinor) = begin
	Spinor(s1.α*s2.α - conj(s2.β)*s1.β,
		   s1.α*s2.β + conj(s2.α)*s1.β)
end

"""
Spinor × Magnetization (Mx + i My, Mz)

V = (x,y,z)

<=>

V = [z -X⃰;

-------- X  -z],	X = x + i y.

Then

(α,β)×(X,z) = (X+,z+),

with

X+ = 2α⃰βz+α⃰²X-β²X⃰

and

z+ = (|α|²-|β|²)z-α⃰ β⃰ X-αβX⃰ .
"""
*(s::Spinor,M::Array) = begin
	[2*conj(s.α)*s.β*M[2]+conj(s.α)^2*M[1]-s.β^2*conj(M[1]),
						   (abs(s.α)^2-abs(s.β)^2)*M[2]-conj(s.α)*conj(s.β)*M[1]-
						   s.α*s.β*conj(M[1])]
end

#Rotation matrices
Rz(φ) = Spinor(exp(-im*φ/2), 0)
Ry(θ) = Spinor(cos(θ/2), sin(θ/2))
"""Rx(θ) = Rz(π/2) Ry(θ) Rz(-π/2)
"""
Rx(θ) = Spinor(cos(θ/2), -im*sin(θ/2))
"""Rg(φ1,θ,φ2) = Rz(φ2) Ry(θ) Rz(φ1)
"""
Rg(φ1,θ,φ2) = Spinor(cos(θ/2)*exp(-im*(φ1+φ2)/2), sin(θ/2)*exp(-im*(φ1-φ2)/2))
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
Q(φ, n::Array{Float64}) = begin
	Spinor(cos(φ/2)-im*n[3]*sin(φ/2),
			-im*(n[1]+im*n[2])*sin(φ/2))
end

"""RF Object"""
struct RF
	A::Complex # Amplitud B1x + i B1y [T]
	T::Float64 # Duration [s]
end
"Duration `T` [s] of the RF array."
dur(x::Array{RF,1}) = sum(x[i].T for i=1:size(x,1))
"Duration `T` [s] of the RF array."
dur(x::Array{RF,2}) = maximum(sum([x[i,j].T for i=1:size(x,1),j=1:size(x,2)],dims=2))

# ## Example
# G = 30e-3
# γ = 2*π*42.5e6
# T = 3e-3; N = 200
# t = 0:T/(N-1):T; Δt = t[2]-t[1]

# B1 = 15e-6; B1e = B1*sinc.(8(t.-T/2)/T)
# φ = -2π*γ*Δt*sqrt.(abs.(B1e).^2) ; φ[φ.==0] .= 1e-17
# n =  2π*γ*Δt*[[B1e[i]/abs(φ[i]) 0 0] for i = 1:N]
# Qs = [Q(φ[i],n[i]) for i=1:N]
# Qt = *(Qs...) #Total rotation matrix

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
# z = .99 * ones(100) * cos.(v)';

# ##
# l1 = PlotlyJS.Layout(;title="Magnetization", yaxis_title="",
#     xaxis_title="t")
# p1 = PlotlyJS.scatter(y=B1e/B1,name="B1e",mode="lines")
# p2 = PlotlyJS.scatter(y=real(Mz),name="Mz",mode="lines")
# p3 = PlotlyJS.scatter(y=imag.(Mxy),name="Mxy",mode="lines")
# pp1 = PlotlyJS.plot([p1,p2,p3],l1)

# l2 = PlotlyJS.Layout(;title="Magnetization",
# 	xaxis_title="Mx",yaxis_title="My",zaxis_tile="Mz")
# p4 = PlotlyJS.scatter3d(x=real(Mxy),y=imag(Mxy),z=real(Mz),xlim=(-1,1),ylim=(-1,1),zlim=(-1,1),
# 	linewidth=5,colorbar=false,legend=false,mode="lines",name="M(t)",line_width=5)
# p5 = PlotlyJS.surface(x=x,y=y,z=z,opacity=.4,showscale=false,colorscale="Greys")
# pp2 = PlotlyJS.plot([p4,p5],l2)

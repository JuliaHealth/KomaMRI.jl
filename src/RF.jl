@doc raw"""
	Spinor(α,β)


Spinor(α,β) with Cayley-Klein parameters α and β.

Based on "Introduction to the Shinnar-Le Roux algorithm", Patrick Le Roux (1995).

A spinor is a way to represent 3D rotations, the underlying representation is a 2 X 2 complex unitary matrix (``\alpha,\beta\in\mathbb{C}``):

```math
R=\left[\begin{array}{cc}
\alpha & -\beta^{*}\\
\beta & \alpha^{*}
\end{array}\right],
```
with ``|\alpha|^2+|\beta|^2 = 1``.

This later operates on the ``2\times2`` representation of ``(x,y,z)`` as follows ``V^{+} = R V R^{*}``.
"""
struct Spinor
	α::Complex
	β::Complex
end
Base.show(io::IO,s::Spinor) = print(io, "Spinor(α = ",round(s.α,digits=3),", β = ",round(s.β,digits=3),")")

"""
Spinor multiplication identity.
	(α1,β1)×(α2,β2) = (α1 α2 - β2⋆ β1 , β2 α1 + α2⋆ β1)
"""
*(s1::Spinor,s2::Spinor) = begin
	Spinor(s1.α*s2.α - conj(s2.β)*s1.β,
		   s1.α*s2.β + conj(s2.α)*s1.β)
end

#Rotation matrices
"Spinor clockwise rotation matrix with angle φ with respect to z-axis."
Rz(φ) = Spinor(exp(-im*φ/2), 0)
"Spinor clockwise rotation matrix with angle Θ with respect to y-axis."
Ry(θ) = Spinor(cos(θ/2), sin(θ/2))
"""
Spinor clockwise rotation matrix with angle Θ with respect to x-axis.

Rx(θ) = Rz(-π/2) Ry(θ) Rz(π/2)
"""
Rx(θ) = Spinor(cos(θ/2), -im*sin(θ/2))
"""
Spinor rotation matrix.

Rg(φ1,θ,φ2) = Rz(φ2) Ry(θ) Rz(φ1)
"""
Rg(φ1,θ,φ2) = Spinor(cos(θ/2)*exp(-im*(φ1+φ2)/2), sin(θ/2)*exp(-im*(φ1-φ2)/2))
"""
Spinor rotation matrix with angle θ with axis in the xy plane u=(cosφ,sinφ).

Rφ(φ,θ) = Rg(-φ,θ,φ) = Rz(φ) Ry(θ) Rz(-φ)
"""
Rφ(φ,θ) = Spinor(cos(θ/2),exp(im*φ)*sin(θ/2))

"""
Pauly, J., Le Roux, P., Nishimura, D., & Macovski, A. (1991).
Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm (NMR imaging).
IEEE Transactions on Medical Imaging, 10(1), 53–65. doi:10.1109/42.75611

Spinor rotation matrix.

Rotation of φ with respect to the axis of rotation n=(nx,ny,nz).

φ = -γ Δt √(|B1|²+(G⋅x)²) = -γ Δt |Bz|

n =  γ Δt/|φ| (B1x, B1y, G⋅x)
"""
Q(φ, n::Array{Float64}) = begin
	Spinor(cos(φ/2)-im*n[3]*sin(φ/2),
		   -im*(n[1]+im*n[2])*sin(φ/2))
end

"""RF Object"""
mutable struct RF
	A::Complex # Amplitud B1x + i B1y [T]
	T::Float64 # Duration [s]
end
"Duration `T` [s] of the RF array Array{RF,1}."
dur(x::Array{RF,1}) = sum(x[i].T for i=1:size(x,1))
"Duration `T` [s] of the RF array Array{RF,2}."
dur(x::Array{RF,2}) = maximum(sum([x[i,j].T for i=1:size(x,1),j=1:size(x,2)],dims=2))
"Generate an RF sequence with amplitudes sampled from a function."
RF_fun(f::Function,T::Real,N::Int64=300) = begin
	RFs = [RF(f(t),T/N) for t = range(0,stop=T,length=N)]
end
getproperty(x::Vector{RF}, f::Symbol) = getproperty.(x,f)
getproperty(x::Matrix{RF}, f::Symbol) = begin
	if f == :A
		x[1,:].A
	elseif f == :T
		x[1,:].T
	end
end
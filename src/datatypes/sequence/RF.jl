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

@doc raw"""
Pauly, J., Le Roux, P., Nishimura, D., & Macovski, A. (1991).
Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm (NMR imaging).
IEEE Transactions on Medical Imaging, 10(1), 53–65. doi:10.1109/42.75611

Spinor rotation matrix.

Rotation of φ with respect to the axis of rotation n=(nx,ny,nz).

```math
\varphi=-\gamma\Delta t\sqrt{\left|B_{1}\right|^{2}+\left(\boldsymbol{G}\cdot\boldsymbol{x}\right)^{2}}=-\gamma\Delta t\left\Vert \boldsymbol{B}\right\Vert
```
```math
\boldsymbol{n}=\boldsymbol{B}/\left\Vert \boldsymbol{B}\right\Vert 
```
"""
Q(φ, nxy, nz) = begin
	Spinor( cos(φ/2)-im*nz*sin(φ/2), -im*nxy*sin(φ/2) )
end
"""
It calculates |\\alpha|^2+|\\beta|^2 of the Cayley-Klein parameters
"""
abs(s::Spinor) = abs(s.α)^2 + abs(s.β)^2


"""RF Object"""
mutable struct RF
	A # Amplitud/Phase B1x + i B1y [T]
	T # Durations [s]
	Δf::Float64    # Frequency offset [Hz]
	delay::Float64 # Delay [s]
	function RF(A,T,Δf,delay)
		@argcheck all(T .>= 0) && delay >= 0 "RF timings must be positive."
		new(A, T, Δf, delay)
    end
	function RF(A,T,Δf)
		@argcheck all(T .>= 0) "RF timings must be positive."
		new(A, T, Δf, 0.)
    end
	function RF(A,T)
		@argcheck all(T .>= 0) >= 0 "RF timings must be positive."
		new(A, T, 0., 0.)
    end
end
#Properties
size(r::RF, i::Int64) = 1 #To fix [r;r;;] concatenation of Julia 1.7.3
*(α::ComplexF64, x::RF) = RF(α*x.A,x.T,x.Δf,x.delay)
"Duration `T` [s] of RF datatype."
dur(x::RF) = sum(x.T)
"Duration `T` [s] of the RF array Array{RF,1}."
dur(x::Array{RF,1}) = sum(sum(x[i].T) for i=1:size(x,1))
"Duration `T` [s] of the RF array Array{RF,2}."
dur(x::Array{RF,2}) = maximum(sum([sum(x[i,j].T) for i=1:size(x,1),j=1:size(x,2)],dims=2))
"Generate an RF sequence with amplitudes sampled from a function."
RF_fun(f::Function,T::Real,N::Int64=300) = begin
	t = range(0,stop=T,length=N)
	A = f.(t)
	RFs = RF(A,T)
end
"Calculates the flip-angle α [deg] of an RF object. α = γ ∫ B1(τ) dτ"
get_flip_angle(x::RF) = begin
	A, NA, T, NT = x.A, length(x.A), x.T, length(x.T)
	dT = T / NA * NT
	α = round(360 * γ * abs(sum(A .* dT)), digits=3) #Pulseq
	α
end
"Calculates the RF center. This calculation includes the RF delay."
get_RF_center(x::RF) = begin
	A, NA, T, NT, delay = x.A, length(x.A), x.T, length(x.T), x.delay
	dT = T / NA * NT .* ones(NA)
	t = cumsum([0; dT])[1:end-1]
	i_center = argmax(abs.(A))
	t_center = t[i_center] + dT[i_center]/2
	t_center + delay
end

one(T::Spinor) = Spinor(1,0)
getproperty(x::Vector{RF}, f::Symbol) = getproperty.(x,f)
getproperty(x::Matrix{RF}, f::Symbol) = begin
	if     f == :Bx
		real.(getproperty.(x,:A))
	elseif f == :By
		imag.(getproperty.(x,:A))
	elseif f == :Δf
		getproperty.(x,:Δf)
	elseif f == :T || f == :delay
		getproperty.(x[1,:],f)
	elseif f == :dur
		T, delay = x.T, x.delay
		ΔT = [sum(t) for t=T] .+ delay
		ΔT
	else
		getproperty.(x,f)
	end
end
#aux
Base.show(io::IO,x::RF) = begin
	r(x) = round.(x,digits=4)
	compact = get(io, :compact, false)
	if !compact
		wave = length(x.A) == 1 ? r(x.A*1e6) : "∿"
		print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"RF($(wave) uT, $(r(sum(x.T)*1e3)) ms, $(r(x.Δf)) Hz)")
	else
		wave = length(x.A) == 1 ? "⊓" : "∿"
		print(io, (sum(abs.(x.A)) > 0 ? wave : "⇿")*"($(r((x.delay+sum(x.T))*1e3)) ms)")
	end
end
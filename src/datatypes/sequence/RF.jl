@doc raw"""
	Spinor(α,β)

Spinor(α,β) with Cayley-Klein parameters α and β. Based on "Introduction to the Shinnar-Le
Roux algorithm", Patrick Le Roux (1995). A spinor is a way to represent 3D rotations, the
underlying representation is a 2 X 2 complex unitary matrix (``\alpha,\beta\in\mathbb{C}``):

```math
R=\left[\begin{array}{cc}
\alpha & -\beta^{*}\\
\beta & \alpha^{*}
\end{array}\right],
```
with ``|\alpha|^2+|\beta|^2 = 1``.

This later operates on the ``2\times2`` representation of ``(x,y,z)`` as follows ``V^{+} =
R V R^{*}``.

# Arguments
- `α::Complex`: the Cayley-Klein parameter α
- `β::Complex`: the Cayley-Klein parameter β
"""
struct Spinor
	α::Complex
	β::Complex
end

"""
    str = show(io::IO,s::Spinor)

Displays the spinnor parameters in the julia REPL.

# Arguments
- `s::Spinor`: the spinnor object

# Returns
- `str` (::String) the output string message
```
"""
Base.show(io::IO,s::Spinor) = begin
    print(io, "Spinor(α = ", round(s.α, digits=3), ", β = ", round(s.β, digits=3), ")")
end

"""
    s = *(s1::Spinor,s2::Spinor)

Spinor multiplication identity: (α1,β1)×(α2,β2) = (α1 α2 - β2⋆ β1 , β2 α1 + α2⋆ β1)

# Arguments
- `s1::Spinor`: first spinnor object
- `s2::Spinor`: second spinnor object

# Returns
- `s::Spinor`: the multiplication spinnor identity result
"""
*(s1::Spinor, s2::Spinor) = begin
	Spinor(s1.α*s2.α - conj(s2.β)*s1.β,
		   s1.α*s2.β + conj(s2.α)*s1.β)
end

"""
    s = Rz(φ)

Spinor clockwise rotation matrix with angle `φ` with respect to z-axis.

# Arguments
- `φ`: (::Real) angle with respect to z-axis

# Returns
- `s::Spinnor`: the spinnor object that represents the `Rz` rotation matrix
"""
Rz(φ) = Spinor(exp(-im*φ/2), 0)

"""
    s = Ry(θ)

Spinor clockwise rotation matrix with angle `θ` with respect to y-axis.

# Arguments
- `θ`: (::Real) angle with respect to y-axis

# Returns
- `s::Spinnor`: the spinnor object that represents the `Ry` rotation matrix
"""
Ry(θ) = Spinor(cos(θ/2), sin(θ/2))

"""
    s = Rx(θ)

Spinor clockwise rotation matrix with angle `θ` with respect to x-axis.

# Arguments
- `θ`: (::Real) angle with respect to x-axis

# Returns
- `s::Spinnor`: the spinnor object that represents the `Rx` rotation matrix
"""
Rx(θ) = Spinor(cos(θ/2), -im*sin(θ/2))

"""
    s = Rg(φ1,θ,φ2)

Spinor rotation matrix: Rg(φ1,θ,φ2) = Rz(φ2) Ry(θ) Rz(φ1)

# Arguments
- `φ1`: (::Real) angle
- `θ`: (::Real) angle
- `φ2`: (::Real) angle

# Returns
- `s::Spinnor`: the spinnor object that represents the `Rg` rotation matrix
"""
Rg(φ1, θ, φ2) = Spinor(cos(θ/2)*exp(-im*(φ1+φ2)/2), sin(θ/2)*exp(-im*(φ1-φ2)/2))

"""
    s = Rφ(φ, θ)

Spinor rotation matrix with angle `θ` with axis in the xy plane u=(cosφ, sinφ).

Rφ(φ,θ) = Rg(-φ,θ,φ) = Rz(φ) Ry(θ) Rz(-φ)

# Arguments
- `φ1`: (::Real) angle
- `θ`: (::Real) angle
- `φ2`: (::Real) angle

# Returns
- `s::Spinnor`: the spinnor object that represents the `Rφ` rotation matrix
"""
Rφ(φ, θ) = Spinor(cos(θ/2), exp(im*φ)*sin(θ/2))

@doc raw"""
    s = Q(φ, nxy, nz)

Spinor rotation matrix. Rotation of `φ` with respect to the axis of rotation n=(nx, ny, nz).

Pauly, J., Le Roux, P., Nishimura, D., & Macovski, A. (1991).
Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm
(NMR imaging).
IEEE Transactions on Medical Imaging, 10(1), 53-65. doi:10.1109/42.75611

```math
\varphi=-\gamma\Delta t\sqrt{\left|B_{1}\right|^{2}+\left(\boldsymbol{G}\cdot\boldsymbol{x}
\right)^{2}}=-\gamma\Delta t\left\Vert \boldsymbol{B}\right\Vert
```
```math
\boldsymbol{n}=\boldsymbol{B}/\left\Vert \boldsymbol{B}\right\Vert
```

# Arguments
- `φ`: (::Real) angle
- `nxy`: (::Real) nxy
- `nz`: (::Real) nz

# Returns
- `s::Spinnor`: the spinnor object that represents the `Q` rotation matrix
"""
Q(φ, nxy, nz) = Spinor(cos(φ/2)-im*nz*sin(φ/2), -im*nxy*sin(φ/2))

"""
    y = abs(s::Spinor)

It calculates |α|^2 + |β|^2 of the Cayley-Klein parameters.

# Arguments
- `s::Spinnor`: the spinnor object

# Returns
- `y`: (::Real) the result of the abs opertor
"""
abs(s::Spinor) = abs(s.α)^2 + abs(s.β)^2


"""
    RF(A, T)
    RF(A, T, Δf)
    RF(A, T, Δf, delay)

The RF Object.

# Arguments
- A: (::Complex{Int64} or ::Vector{Complex{Int64}}) Amplitud/Phase B1x + i B1y [T]
- T: (::Int64 or ::Vector{Int64}) duurations in [s]
- `Δf::Float64`: the frequency offset in [Hz]
- `delay::Float64`: the delay time in [s]
"""
mutable struct RF
	A                # Amplitud/Phase B1x + i B1y [T]
	T                # Durations [s]
	Δf::Float64      # Frequency offset [Hz]
	delay::Float64   # Delay [s]
	function RF(A, T, Δf, delay)
		@argcheck all(T .>= 0) && delay >= 0 "RF timings must be positive."
		new(A, T, Δf, delay)
    end
	function RF(A, T, Δf)
		@argcheck all(T .>= 0) "RF timings must be positive."
		new(A, T, Δf, 0.)
    end
	function RF(A, T)
		@argcheck all(T .>= 0) >= 0 "RF timings must be positive."
		new(A, T, 0., 0.)
    end
end

"""
    str = show(io::IO, x::RF)

Displays information about the RF object `x` in the julia REPL.

# Arguments
- `x::RF`: the RF object

# Returns
- `str` (::String) the output string message
"""
Base.show(io::IO, x::RF) = begin
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

"""
    y = getproperty(x::Vector{RF}, f::Symbol)
    y = getproperty(x::Matrix{RF}, f::Symbol)

Overchages Base.getproperty(). It is meant to access properties of the RF vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x::Vector{RF}`: the vector of RF objects
- `x::Matrix{RF}`: the matrix of RF objects
- `f::Symbol`: custom options are the `:Bx`, `:By`, `:Δf`, `:T`, `:delay` and `:dur` symbols

# Returns
- `y`: (::Vector{Any} or ::Matrix{Any}) the vector with the property defined by the symbol
    `f` for all elements of the RF vector or matrix `x`
"""
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

# Properties
one(T::Spinor) = Spinor(1,0)
size(r::RF, i::Int64) = 1 #To fix [r;r;;] concatenation of Julia 1.7.3
*(α::ComplexF64, x::RF) = RF(α*x.A,x.T,x.Δf,x.delay)

"""
   y = dur(x::RF)
   y = dur(x::Array{RF,1})
   y = dur(x::Array{RF,2})

Duration time in [s] of RF object or RF array.

# Arguments
- `x`: (::RF, ::Array{RF,1} or ::Array{RF,2}) the RF object or RF array

# Returns
- `y`: (::Float64) the duration of the RF object or RF array in [s]
"""
dur(x::RF) = sum(x.T)
dur(x::Array{RF,1}) = sum(sum(x[i].T) for i=1:size(x,1))
dur(x::Array{RF,2}) = maximum(sum([sum(x[i,j].T) for i=1:size(x,1),j=1:size(x,2)],dims=2))

"""
    rf = RF_fun(f::Function, T::Real, N::Int64)

Generate an RF sequence with amplitudes sampled from a function.

!!! note
    This function is not being used in this KomaMRI version.

# Arguments
- `f::Function`: the function for the RF amplitud
- `T::Real`: the duration of the RF pulse
- `N::Int64`: the number of samples of the RF pulse

# Returns
- `rf::RF`: the rf object with aplitud defined by the function `f`
"""
RF_fun(f::Function, T::Real, N::Int64=300) = begin
	t = range(0, stop=T, length=N)
	A = f.(t)
	RF(A, T)
end

"""
    α = get_flip_angle(x::RF)

Calculates the flip angle α [deg] of an RF object. α = γ ∫ B1(τ) dτ

# Arguments
- `x::RF`: the RF pulse

# Returns
- `α`: (::Int64) the flip angle RF pulse `x` in [deg]
"""
get_flip_angle(x::RF) = begin
	A, NA, T, NT = x.A, length(x.A), x.T, length(x.T)
	dT = T / NA * NT
	α = round(360 * γ * abs(sum(A .* dT)), digits=3) #Pulseq
	return α
end

"""
    t = get_RF_center(x::RF)

Calculates the time where is the center of the RF pulse `x`. This calculation includes the
RF delay.

# Arguments
- `x::RF`: the RF pulse

# Returns
- `t`: (::Int64) the time where is the center of the RF pulse `x`
"""
get_RF_center(x::RF) = begin
	A, NA, T, NT, delay = x.A, length(x.A), x.T, length(x.T), x.delay
	dT = T / NA * NT .* ones(NA)
	t = cumsum([0; dT])[1:end-1]
	i_center = argmax(abs.(A))
	t_center = t[i_center] + dT[i_center]/2
	return t_center + delay
end

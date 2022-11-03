@doc raw"""
	spinor = Spinor(α, β)

Spinor(α, β) with Cayley-Klein parameters α and β. Based on "Introduction to the Shinnar-Le
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
- `α`: (`::Complex{Int64}`) the Cayley-Klein parameter α
- `β`: (`::Complex{Int64}`) the Cayley-Klein parameter β

# Returns
- `spinor`: (`::Spinor`) the Spinor struct
"""
struct Spinor
	α::Complex
	β::Complex
end

"""
    str = show(io::IO, s::Spinor)

Displays the spinor parameters in the julia REPL.

# Arguments
- `s`: (`::Spinor`) the spinnor struct

# Returns
- `str`: (`::String`) the output string message
"""
Base.show(io::IO, s::Spinor) = begin
    print(io, "Spinor(α = ", round(s.α, digits=3), ", β = ", round(s.β, digits=3), ")")
end

"""
    s = *(s1::Spinor, s2::Spinor)

Spinor multiplication identity: (α1,β1)×(α2,β2) = (α1 α2 - β2⋆ β1 , β2 α1 + α2⋆ β1)

# Arguments
- `s1`: (`::Spinor`) the first spinnor struct
- `s2`: (`::Spinor`) the second spinnor struct

# Returns
- `s`: (`::Spinor`) the multiplication spinnor identity result
"""
*(s1::Spinor, s2::Spinor) = begin
	Spinor(s1.α*s2.α - conj(s2.β)*s1.β,
		   s1.α*s2.β + conj(s2.α)*s1.β)
end

"""
    s = Rz(φ)

Spinor clockwise rotation matrix with angle `φ` with respect to z-axis.

# Arguments
- `φ`: (`::Real`, `[rad]`) the angle with respect to z-axis

# Returns
- `s`: (`::Spinnor`) the spinnor struct that represents the `Rz` rotation matrix
"""
Rz(φ) = Spinor(exp(-im*φ/2), 0)

"""
    s = Ry(θ)

Spinor clockwise rotation matrix with angle `θ` with respect to y-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) the angle with respect to y-axis

# Returns
- `s`: (`::Spinnor`) the spinnor struct that represents the `Ry` rotation matrix
"""
Ry(θ) = Spinor(cos(θ/2), sin(θ/2))

"""
    s = Rx(θ)

Spinor clockwise rotation matrix with angle `θ` with respect to x-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) the angle with respect to x-axis

# Returns
- `s`: (`::Spinnor`) the spinnor struct that represents the `Rx` rotation matrix
"""
Rx(θ) = Spinor(cos(θ/2), -im*sin(θ/2))

"""
    s = Rg(φ1, θ, φ2)

Spinor rotation matrix: Rg(φ1, θ, φ2) = Rz(φ2) Ry(θ) Rz(φ1)

# Arguments
- `φ1`: (`::Real`, `[rad]`) the φ1 angle
- `θ`: (`::Real`, `[rad]`) the θ angle
- `φ2`: (`::Real`, `[rad]`) the φ2 angle

# Returns
- `s`: (`::Spinnor`) the spinnor struct that represents the `Rg` rotation matrix
"""
Rg(φ1, θ, φ2) = Spinor(cos(θ/2)*exp(-im*(φ1+φ2)/2), sin(θ/2)*exp(-im*(φ1-φ2)/2))

"""
    s = Rφ(φ, θ)

Spinor rotation matrix with angle `θ` with axis in the xy plane u=(cosφ, sinφ).

Rφ(φ,θ) = Rg(-φ,θ,φ) = Rz(φ) Ry(θ) Rz(-φ)

# Arguments
- `φ`: (`::Real`, `[rad]`) the φ angle
- `θ`: (`::Real`, `[rad]`) the θ angle

# Returns
- `s`: (`::Spinnor`) the spinnor struct that represents the `Rφ` rotation matrix
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
- `φ`: (`::Real`, `[rad]`) the φ angle
- `nxy`: (`::Real`) the nxy factor
- `nz`: (`::Real`) the nz factor

# Returns
- `s`: (`::Spinnor`) the spinnor struct that represents the `Q` rotation matrix
"""
Q(φ, nxy, nz) = Spinor(cos(φ/2)-im*nz*sin(φ/2), -im*nxy*sin(φ/2))

"""
    y = abs(s::Spinor)

It calculates |α|^2 + |β|^2 of the Cayley-Klein parameters.

# Arguments
- `s`: (`::Spinnor`) the spinnor struct

# Returns
- `y`: (`::Real`) the result of the abs operator
"""
abs(s::Spinor) = abs(s.α)^2 + abs(s.β)^2


"""
    rf = RF(A, T)
    rf = RF(A, T, Δf)
    rf = RF(A, T, Δf, delay)

The RF struct.

# Arguments
- `A`: (`::Complex{Int64}`, `[T]`) the amplitud-phase B1x + i B1y
- `T`: (`::Int64`, [`s`]) the durations of the RF
- `Δf`: (`::Float64`, [`Hz`]) the frequency offset of the RF
- `delay`: (`::Float64`, [`s`]) the delay time of the RF

# Returns
- `rf`: (`::RF`) the RF struct

# Examples
```julia-repl
julia> d1, d2, d3 = 0.8, 0.4, 0.8;

julia> fsinc = x -> 2 * sinc(3*pi*(x - d1/2)) * 1e-3;

julia> matrixGrads = [Grad(0, d1) Grad(0, d2) Grad(0, d3)];

julia> matrixRFs = [KomaMRI.RF_fun(fsinc, d1) RF(0, d2) RF(0, d3)];

julia> seq = Sequence(matrixGrads, matrixRFs)
Sequence[ τ = 2000.0 ms | blocks: 3 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(seq)
```
"""
mutable struct RF
	A
	T
	Δf::Float64
	delay::Float64
	function RF(A, T, Δf, delay)
        any(T .< 0) || delay < 0 ? error("RF timings must be non-negative.") : new(A, T, Δf, delay)
    end
	function RF(A, T, Δf)
        any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, Δf, 0.)
    end
	function RF(A, T)
        any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, 0., 0.)
    end
end

"""
    str = show(io::IO, x::RF)

Displays information about the RF struct `x` in the julia REPL.

# Arguments
- `x`: (`::RF`) the RF struct

# Returns
- `str`: (`::String`) the output string message
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
- `x`: (`::Vector{RF}` or `::Matrix{RF}`) the vector or matrix of RF structs
- `f`: (`::Symbol`, opts: [`:A`, `:Bx`, `:By`, `:T`, `:Δf`, `:delay` and `:dur`]) the input
    symbol that represents a property of the vector or matrix of RF structs

# Returns
- `y`: (`::Vector{Any}` or `::Matrix{Any}`) the vector with the property defined by the
    symbol `f` for all elements of the RF vector or matrix `x`
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

Duration time in [s] of RF struct or RF array.

# Arguments
- `x`: (`::RF` or `::Array{RF,1}` or `::Array{RF,2}`) the RF struct or RF array

# Returns
- `y`: (`::Float64`, [`s`]) the duration of the RF struct or RF array
"""
dur(x::RF) = sum(x.T)
dur(x::Array{RF,1}) = sum(sum(x[i].T) for i=1:size(x,1))
dur(x::Array{RF,2}) = maximum(sum([sum(x[i,j].T) for i=1:size(x,1),j=1:size(x,2)],dims=2))

"""
    rf = RF_fun(f::Function, T::Real, N::Int64)

Generate an RF sequence with amplitudes sampled from a function waveform.

!!! note
    This function is not being used in this KomaMRI version.

# Arguments
- `f`: (`::Function`, [`T`]) the function for the RF amplitud waveform
- `T`: (`::Real`, [`s`]) the duration of the RF pulse
- `N`: (`::Int64`) the number of samples of the RF pulse

# Returns
- `rf`:(`::RF`) the RF struct with amplitud defined by the function `f`
"""
RF_fun(f::Function, T::Real, N::Int64=300) = begin
	t = range(0, stop=T, length=N)
	A = f.(t)
	RF(A, T)
end

"""
    α = get_flip_angle(x::RF)

Calculates the flip angle α [deg] of an RF struct. α = γ ∫ B1(τ) dτ

# Arguments
- `x`: (`::RF`) the RF struct

# Returns
- `α`: (`::Int64`, `[deg]`) the flip angle RF struct `x`
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
- `x`: (`::RF`) the RF struct

# Returns
- `t`: (`::Int64`, `[s]`) the time where is the center of the RF pulse `x`
"""
get_RF_center(x::RF) = begin
	A, NA, T, NT, delay = x.A, length(x.A), x.T, length(x.T), x.delay
	dT = T / NA * NT .* ones(NA)
	t = cumsum([0; dT])[1:end-1]
	i_center = argmax(abs.(A))
	t_center = t[i_center] + dT[i_center]/2
	return t_center + delay
end

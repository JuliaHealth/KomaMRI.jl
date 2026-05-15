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
- `α`: (`::Complex{Float64}`) Cayley-Klein parameter α
- `β`: (`::Complex{Float64}`) Cayley-Klein parameter β

# Returns
- `spinor`: (`::Spinor`) Spinor struct
"""
struct Spinor{T<:Real,V<:AbstractVector}
	α::V
	β::V
end
function Spinor(α::V, β::V) where {V<:AbstractVector}
	T = promote_type(_real_storage_eltype(typeof(α)), _real_storage_eltype(typeof(β)))
	return Spinor{T,V}(α, β)
end
Spinor(α::Complex{T}, β::Complex{T}) where {T<:Real} = Spinor([α], [β])
Spinor(α::T, β::T) where {T<:Real} = Spinor([complex(α)], [complex(β)])
one(T::Spinor) = Spinor(1.,0.)
Base.getindex(s::Spinor, i) = Spinor(s.α[i], s.β[i])
Base.view(s::Spinor, i::UnitRange) = @views Spinor(s.α[i], s.β[i])

function _storage_eltype(::Type{A}) where {A}
    params = Base.unwrap_unionall(A).parameters
    return !isempty(params) && params[1] isa Type ? params[1] : eltype(A)
end
_real_storage_eltype(::Type{A}) where {A} = _real_eltype(_storage_eltype(A))
_real_eltype(::Type{Complex{T}}) where {T} = T
_real_eltype(::Type{T}) where {T<:Real} = T
"""
    str = show(io::IO, s::Spinor)

Displays the spinor parameters in the julia REPL.

# Arguments
- `s`: (`::Spinor`) Spinor struct

# Returns
- `str`: (`::String`) output string message
"""
Base.show(io::IO, s::Spinor) = begin
    print(io, "Spinor(α = ", s.α, ", β = ", s.β, ")")
end

"""
    s = *(s2::Spinor, s1::Spinor)

Spinor multiplication identity: (α2,β2) × (α1,β1) = (α1 α2 - β2⋆ β1 , β2 α1 + α2⋆ β1)
"""
*(s2::Spinor, s1::Spinor) = begin
	Spinor(s1.α.*s2.α .- conj.(s2.β).*s1.β,
		   s1.α.*s2.β .+ conj.(s2.α).*s1.β)
end

"""
    s = Rz(φ)

Spinor counter-clockwise rotation matrix with angle `φ` with respect to z-axis.

# Arguments
- `φ`: (`::Real`, `[rad]`) angle with respect to z-axis

# Returns
- `s`: (`::Spinor`) spinnor struct that represents the `Rz` rotation matrix
"""
Rz(φ) = Spinor(exp(-1im*φ/2), 0.0im)

"""
    s = Ry(θ)

Spinor counter-clockwise rotation matrix with angle `θ` with respect to y-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) angle with respect to y-axis

# Returns
- `s`: (`::Spinor`) spinor struct that represents the `Ry` rotation matrix
"""
Ry(θ) = Spinor(cos(θ/2)+0im, sin(θ/2)+0im)

"""
    s = Rx(θ)

Spinor counter-clockwise rotation matrix with angle `θ` with respect to x-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) angle with respect to x-axis

# Returns
- `s`: (`::Spinor`) spinor struct that represents the `Rx` rotation matrix
"""
Rx(θ) = Spinor(cos(θ/2)+0im, -1im*sin(θ/2))

"""
    s = Rφ(φ, θ)

Spinor counter-clockwise rotation matrix of angle `θ` around an axis in the xy plane.
The rotation axis makes an angle `φ` with the axis y, or u=(sinφ, cosφ, 0).

Rφ(φ,θ) = Rg(-φ,θ,φ) = Rz(-φ) Ry(θ) Rz(φ)

# Arguments
- `φ`: (`::Real`, `[rad]`) φ angle
- `θ`: (`::Real`, `[rad]`) θ angle

# Returns
- `s`: (`::Spinor`) spinnor struct that represents the `Rφ` rotation matrix
"""
Rφ(φ, θ) = Spinor(cos(θ/2)+0.0im, exp(-1im*φ)*sin(θ/2))

"""
    s = Rg(φ1, θ, φ2)

General Spinor rotation matrix: Rg(φ1, θ, φ2) = Rz(φ2) Ry(θ) Rz(φ1)

# Arguments
- `φ1`: (`::Real`, `[rad]`) φ1 angle
- `θ`: (`::Real`, `[rad]`) θ angle
- `φ2`: (`::Real`, `[rad]`) φ2 angle

# Returns
- `s`: (`::Spinor`) spinor struct that represents the `Rg` rotation matrix
"""
Rg(φ1, θ, φ2) = Spinor(cos(θ/2)*exp(-1im*(φ1+φ2)/2), sin(θ/2)*exp(-1im*(φ1-φ2)/2))


@doc raw"""
    s = Q(φ, nxy, nz)

Spinor rotation matrix. Counter-clockwise rotation of `φ` with respect to the axis of rotation n=(nx, ny, nz).

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
- `φ`: (`::Real`, `[rad]`) φ angle
- `nxy`: (`::Real`) nxy factor
- `nz`: (`::Real`) nz factor

# Returns
- `s`: (`::Spinor`) spinnor struct that represents the `Q` rotation matrix
"""
function Q(φ, nxy, nz)
    φ_half = φ ./ 2
    sin_φ_half = sin.(φ_half)
    neg_im = complex.(zero.(φ), -Base.one.(φ))
    return Spinor(
        cos.(φ_half) .+ neg_im .* nz .* sin_φ_half,
        neg_im .* nxy .* sin_φ_half,
    )
end

"""
    y = abs(s::Spinor)

It calculates |α|^2 + |β|^2 of the Cayley-Klein parameters.

# Arguments
- `s`: (`::Spinor`) spinnor struct

# Returns
- `y`: (`::Real`) result of the abs operator
"""
abs(s::Spinor) = abs.(s.α).^2 + abs.(s.β).^2

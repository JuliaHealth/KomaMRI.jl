@doc raw"""
    translate = Translate(dx, dy, dz)

Translate struct. It produces a linear translation.
Its fields are the final displacements in the three axes (dx, dy, dz).

# Arguments
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z

# Returns
- `translate`: (`::Translate`) Translate struct

# Examples
```julia-repl
julia> translate = Translate(dx=0.01, dy=0.02, dz=0.03)
```
"""
@with_kw struct Translate{T<:Real} <: SimpleAction{T}
    dx         :: T
    dy         :: T
    dz         :: T
end

TranslateX(dx::T) where {T<:Real} = Translate(dx, zero(T), zero(T))
TranslateY(dy::T) where {T<:Real} = Translate(zero(T), dy, zero(T))
TranslateZ(dz::T) where {T<:Real} = Translate(zero(T), zero(T), dz)

function displacement_x!(
    ux::AbstractArray{T},
    action::Translate{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    ux .= t.* action.dx
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
    action::Translate{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    uy .= t .* action.dy
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    action::Translate{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    uz .= t .* action.dz
    return nothing
end


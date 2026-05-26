@doc raw"""
    t = Translate(dx, dy, dz)

Translate struct. It produces a linear translation.
Its fields are the final displacements in the three axes (dx, dy, dz).

# Arguments
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z

# Returns
- `t`: (`::Translate`) Translate struct

# Examples
```julia-repl
julia> t = Translate(dx=0.01, dy=0.02, dz=0.03)
```
"""
@with_kw struct Translate <: SimpleAction
    dx :: Real
    dy :: Real
    dz :: Real
end

TranslateX(dx::Real) = Translate(dx, zero(dx), zero(dx))
TranslateY(dy::Real) = Translate(zero(dy), dy, zero(dy))
TranslateZ(dz::Real) = Translate(zero(dz), zero(dz), dz)

is_composable(m::Translate) = false

function displacement_x!(ux, action::Translate, x, y, z, t)
    ux .= t.* action.dx
    return nothing
end

function displacement_y!(uy, action::Translate, x, y, z, t)
    uy .= t .* action.dy
    return nothing
end

function displacement_z!(uz, action::Translate, x, y, z, t)
    uz .= t .* action.dz
    return nothing
end


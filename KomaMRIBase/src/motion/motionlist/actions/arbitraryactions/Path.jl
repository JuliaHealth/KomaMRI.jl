@doc raw"""
    path = Path(dx, dy, dz)

Path motion struct. (...)

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z

# Returns
- `path`: (`::Path`) Path struct

# Examples
```julia-repl
julia> p = Path(dx=[0.01 0.02], dy=[0.02 0.03], dz=[0.03 0.04])
```
"""
@with_kw struct Path{T<:Real} <: ArbitraryAction{T}
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
end
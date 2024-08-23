@doc raw"""
    flowpath = FlowPath(dx, dy, dz)

FlowPath motion struct. (...)

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `spin_reset`: (`::AbstractArray{Bool}`) reset spin state flags

# Returns
- `flowpath: (`::FlowPath`) FlowPath struct

# Examples
```julia-repl
julia> fp = FlowPath(dx=[0.01 0.02], dy=[0.02 0.03], dz=[0.03 0.04], spin_reset=[false, false])
```
"""
@with_kw struct FlowPath{T<:Real} <: ArbitraryAction{T}
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
    spin_reset::AbstractArray{Bool}
end
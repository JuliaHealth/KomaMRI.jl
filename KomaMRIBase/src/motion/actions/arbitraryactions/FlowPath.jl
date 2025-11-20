@doc raw"""
    f = FlowPath(dx, dy, dz, spin_reset)

FlowPath struct. This action is the same as `Path`, 
except that it includes an additional field, called `spin_reset`, 
which accounts for spins leaving the volume and being remapped 
to another input position. When this happens, the magnetization 
state of these spins must be reset during the simulation. 

As with the `dx`, `dy` and `dz` matrices, `spin_reset`
has a size of (``N_{spins} \times \; N_{discrete\,times}``).

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `spin_reset`: (`::AbstractArray{Bool}`) reset spin state flags

# Returns
- `f`: (`::FlowPath`) FlowPath struct

# Examples
```julia-repl
julia> f = FlowPath(
           dx=[0.01 0.02; 0.02 0.03], 
           dy=[0.02 0.03; 0.03 0.04], 
           dz=[0.03 0.04; 0.04 -0.04],
           spin_reset=[false false; false true]
       )
```
"""
@with_kw struct FlowPath{T<:Real} <: ArbitraryAction{T}
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
    spin_reset::AbstractArray{Bool}
end

FlowPath(dx::AbstractArray{T}, dy::AbstractArray{T}, dz::AbstractArray{T}, spin_reset::BitMatrix) where T<:Real = FlowPath(dx, dy, dz, collect(spin_reset))

function add_reset_times!(t, a::FlowPath, t_start, t_end, periods)
    aux = t_start .+ (t_end - t_start)/(size(a.spin_reset)[2]-1) * (getindex.(findall(a.spin_reset .== 1), 2) .- 1)
    append!(t, times(aux, t_start, t_end, periods) .- MIN_RISE_TIME)
end
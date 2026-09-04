@doc raw"""
    f = FlowPath(dx, dy, dz, spin_reset; cycle_map=nothing)

FlowPath struct. This action is the same as `Path`, 
except that it includes an additional field, called `spin_reset`, 
which accounts for spins leaving the volume and being remapped 
to another input position. When this happens, the magnetization 
state of these spins must be reset during the simulation. 

As with the `dx`, `dy` and `dz` matrices, `spin_reset`
has a size of (``N_{spins} \times \; N_{discrete\,times}``).

For periodic, non-closed trajectories, `cycle_map` can provide a precomputed
destination-to-source particle mapping. At each cycle boundary, particle `i`
receives the magnetization of particle `cycle_map[i]`. The map uses 1-based
indices local to this `FlowPath`. Koma does not calculate this map. A map is only
meaningful for the whole particle set, so slicing a `FlowPath` drops it.

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `spin_reset`: (`::AbstractArray{Bool}`) reset spin state flags
- `cycle_map`: (`::Union{Nothing,AbstractVector{Int}}`) optional periodic magnetization map

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
    cycle_map::Union{Nothing,Vector{Int}} = nothing
end

FlowPath(dx, dy, dz, spin_reset::BitMatrix; cycle_map=nothing) = FlowPath(dx, dy, dz, collect(spin_reset); cycle_map)
function FlowPath(dx, dy, dz, spin_reset::AbstractArray{Bool}; cycle_map=nothing)
    isnothing(cycle_map) ||
        (eltype(cycle_map) !== Bool && length(cycle_map) == size(dx, 1) && all(in(axes(dx, 1)), cycle_map)) ||
        throw(ArgumentError("cycle_map must hold one source particle index per FlowPath particle."))
    return FlowPath(dx, dy, dz, spin_reset, isnothing(cycle_map) ? nothing : Vector{Int}(cycle_map))
end

# A cycle_map indexes the whole particle set, so a sub-group cannot carry it.
_sliced_cycle_map(a::FlowPath, p) = (p isa Colon || p == 1:size(a.dx, 1)) ? a.cycle_map : nothing

Base.getindex(a::FlowPath, p) = FlowPath(a.dx[p, :], a.dy[p, :], a.dz[p, :], a.spin_reset[p, :], _sliced_cycle_map(a, p))
Base.view(a::FlowPath, p) = @views FlowPath(a.dx[p, :], a.dy[p, :], a.dz[p, :], a.spin_reset[p, :], _sliced_cycle_map(a, p))

function add_reset_times!(t, a::FlowPath, t_start, t_end, periods)
    aux = t_start .+ (t_end - t_start)/(size(a.spin_reset)[2]-1) * (getindex.(findall(a.spin_reset .== 1), 2) .- 1)
    append!(t, times(aux, t_start, t_end, periods) .- MIN_RISE_TIME)
end

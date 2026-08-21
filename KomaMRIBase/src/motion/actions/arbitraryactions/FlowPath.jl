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
indices local to this `FlowPath`. Koma does not calculate this map. Remapping is
currently supported only on the GPU, and phantoms using it cannot be sliced.

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
@with_kw struct FlowPath{T<:Real,M} <: ArbitraryAction{T}
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
    spin_reset::AbstractArray{Bool}
    cycle_map::M = nothing
end

function FlowPath(
    dx::AbstractArray{T},
    dy::AbstractArray{T},
    dz::AbstractArray{T},
    spin_reset::AbstractArray{Bool},
    cycle_map::Union{Nothing,AbstractVector{<:Integer}},
) where T<:Real
    size(dx) == size(dy) == size(dz) == size(spin_reset) ||
        throw(DimensionMismatch("FlowPath arrays must have equal sizes."))
    if !isnothing(cycle_map)
        eltype(cycle_map) === Bool && throw(ArgumentError("cycle_map must contain particle indices, not Bool values."))
        length(cycle_map) == size(dx, 1) ||
            throw(DimensionMismatch("cycle_map length must equal the number of FlowPath particles."))
        all(i -> 1 <= i <= size(dx, 1), cycle_map) ||
            throw(ArgumentError("cycle_map indices must be between 1 and $(size(dx, 1))."))
        cycle_map = Int.(collect(cycle_map))
    end
    spin_reset = spin_reset isa BitMatrix ? collect(spin_reset) : spin_reset
    return FlowPath{T,typeof(cycle_map)}(dx, dy, dz, spin_reset, cycle_map)
end

function FlowPath(
    dx::AbstractArray{T},
    dy::AbstractArray{T},
    dz::AbstractArray{T},
    spin_reset::AbstractArray{Bool};
    cycle_map=nothing,
) where T<:Real
    return FlowPath(dx, dy, dz, spin_reset, cycle_map)
end

function Base.getindex(action::FlowPath{T}, p) where T
    isnothing(action.cycle_map) || _is_full_spin_selection(p, size(action.dx, 1)) ||
        throw(ArgumentError("FlowPath with cycle_map cannot be sliced."))
    return FlowPath(
        action.dx[p, :], action.dy[p, :], action.dz[p, :], action.spin_reset[p, :],
        action.cycle_map,
    )
end

function Base.view(action::FlowPath{T}, p) where T
    isnothing(action.cycle_map) || _is_full_spin_selection(p, size(action.dx, 1)) ||
        throw(ArgumentError("FlowPath with cycle_map cannot be sliced."))
    return @views FlowPath{T,typeof(action.cycle_map)}(
        action.dx[p, :], action.dy[p, :], action.dz[p, :], action.spin_reset[p, :],
        action.cycle_map,
    )
end

function add_reset_times!(t, a::FlowPath, t_start, t_end, periods)
    aux = t_start .+ (t_end - t_start)/(size(a.spin_reset)[2]-1) * (getindex.(findall(a.spin_reset .== 1), 2) .- 1)
    append!(t, times(aux, t_start, t_end, periods) .- MIN_RISE_TIME)
end

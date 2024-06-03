# ------ SimpleMotionType
abstract type SimpleMotionType{T<:Real} end

is_composable(motion_type::SimpleMotionType{T}) where {T<:Real} = false

"""
    motion = SimpleMotion(types)

SimpleMotion model. It allows for the definition of motion by means of simple parameters.
The `SimpleMotion` struct is composed by only one field, called `types`, 
which is a vector of simple motion types. This vector will contain as many elements
as simple motions we want to combine.

# Arguments
- `types`: (`::Vector{<:SimpleMotionType{T}}`) vector of simple motion types

# Returns
- `motion`: (`::SimpleMotion`) SimpleMotion struct

# Examples
```julia-repl
julia> motion = SimpleMotion([
            Translation(dx=0.01, dy=0.02, dz=0.0, t_start=0.0, t_end=0.5),
            Rotation(pitch=15.0, roll=0.0, yaw=20.0, t_start=0.1, t_end=0.5),
            HeartBeat(circumferential_strain=-0.3, radial_strain=-0.2, longitudinal_strain=0.0, t_start=0.2, t_end=0.5)
        ])
```
"""
struct SimpleMotion{T<:Real} <: MotionModel{T}
    types::Vector{<:SimpleMotionType{T}}
end

Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
function Base.getindex(
    motion::SimpleMotion,
    p::Union{AbstractRange,AbstractVector,Colon},
    q::Union{AbstractRange,AbstractVector,Colon},
)
    return motion
end

Base.vcat(m1::SimpleMotion, m2::SimpleMotion) = SimpleMotion([m1.types; m2.types])

Base.:(==)(m1::SimpleMotion, m2::SimpleMotion) = reduce(&, [m1.types[i] == m2.types[i] for i in 1:length(m1.types)])
Base.:(≈)(m1::SimpleMotion, m2::SimpleMotion)  = reduce(&, [m1.types[i] ≈ m2.types[i] for i in 1:length(m1.types)])
# When they are the same type  
Base.:(==)(t1::T, t2::T) where {T<:SimpleMotionType} = reduce(&, [getfield(t1, field) == getfield(t2, field) for field in fieldnames(T)])
Base.:(≈)(t1::T, t2::T) where {T<:SimpleMotionType}  = reduce(&, [getfield(t1, field) ≈ getfield(t2, field) for field in fieldnames(T)])
# When they are not (default)  
Base.:(==)(t1::SimpleMotionType, t2::SimpleMotionType) = false
Base.:(≈)(t1::SimpleMotionType, t2::SimpleMotionType)  = false

"""
    x, y, z = get_spin_coords(motion, x, y, z, t)

Calculates the position of each spin at a set of arbitrary time instants, i.e. the time steps of the simulation. 
For each dimension (x, y, z), the output matrix has ``N_{\text{spins}}`` rows and `length(t)` columns.

# Arguments
- `motion`: (`::MotionModel`) phantom motion
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `t`: (`::AbstractArray{T<:Real}`) horizontal array of time instants

# Returns
- `x, y, z`: (`::Tuple{AbstractArray, AbstractArray, AbstractArray}`) spin positions over time
"""
function get_spin_coords(
    motion::SimpleMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    xt, yt, zt = x .+ 0*t, y .+ 0*t, z .+ 0*t
    # Composable motions: they need to be run sequentially
    for motion in Iterators.filter(is_composable, motion.types)
        aux = xt .+ 0, yt .+ 0, zt .+ 0
        xt .+= displacement_x(motion, aux..., t)
        yt .+= displacement_y(motion, aux..., t)
        zt .+= displacement_z(motion, aux..., t)
    end
    # Additive motions: these motions can be run in parallel
    for motion in Iterators.filter(!is_composable, motion.types)
        xt .+= displacement_x(motion, x, y, z, t)
        yt .+= displacement_y(motion, x, y, z, t)
        zt .+= displacement_z(motion, x, y, z, t)
    end
    return xt, yt, zt
end

function times(motion::SimpleMotion)
    nodes = reduce(vcat, [times(type) for type in motion.types])
    nodes = unique(sort(nodes))
    return nodes
end

function sort_motions!(motion::SimpleMotion)
    return sort!(motion.types; by=m -> times(m)[1])
end

# --------- Simple Motion Types: -------------
# Non-periodic types: defined by an initial time (t_start), an end time (t_end) and a displacement      
include("simplemotion/Translation.jl")
include("simplemotion/Rotation.jl")
include("simplemotion/HeartBeat.jl")

function unit_time(t::AbstractArray{T}, t_start::T, t_end::T) where {T<:Real}
    if t_start == t_end
        return (t .>= t_start) .* one(T) # The problem with this is that it returns a BitVector (type stability issues)
    else
        return min.(max.((t .- t_start) ./ (t_end - t_start), zero(T)), one(T))
    end
end

# Periodic types: defined by the period, the temporal symmetry and a displacement (amplitude)
include("simplemotion/PeriodicTranslation.jl")
include("simplemotion/PeriodicRotation.jl")
include("simplemotion/PeriodicHeartBeat.jl")

function unit_time_triangular(t, period, asymmetry)
    t_rise = period * asymmetry
    t_fall = period * (1 - asymmetry)
    t_relative = mod.(t, period)
    t_unit =
        ifelse.(
            t_relative .< t_rise,
            t_relative ./ t_rise,
            1 .- (t_relative .- t_rise) ./ t_fall,
        )
    return t_unit
end

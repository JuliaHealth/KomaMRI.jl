# ------ SimpleMotionType
abstract type SimpleMotionType{T<:Real} end

is_composable(motion_type::SimpleMotionType{T}) where {T<:Real} = false

"""
    motion = SimpleMotion(types)

SimpleMotion model. It allows for the definition of motion by means of simple parameters.
The `SimpleMotion` struct is composed by only one field, called `types`, 
which is a tuple of simple motion types. This tuple will contain as many elements
as simple motions we want to combine.

# Arguments
- `types`: (`::Tuple{Vararg{<:SimpleMotionType{T}}}`) vector of simple motion types

# Returns
- `motion`: (`::SimpleMotion`) SimpleMotion struct

# Examples
```julia-repl
julia> motion = SimpleMotion(
            Translation(dx=0.01, dy=0.02, dz=0.0, t_start=0.0, t_end=0.5),
            Rotation(pitch=15.0, roll=0.0, yaw=20.0, t_start=0.1, t_end=0.5),
            HeartBeat(circumferential_strain=-0.3, radial_strain=-0.2, longitudinal_strain=0.0, t_start=0.2, t_end=0.5)
        )
```
"""
struct SimpleMotion{T<:Real} <: MotionModel{T}
    types::Tuple{Vararg{<:SimpleMotionType{T}}}
end

SimpleMotion(types...) = SimpleMotion(types)

Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.view(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion

Base.vcat(m1::SimpleMotion, m2::SimpleMotion) = SimpleMotion(m1.types..., m2.types...)

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
    t::AbstractArray{T}
) where {T<:Real}
    motion = sort_motion(motion)
    # Buffers for positions:
    xt, yt, zt = x .+ 0*t, y .+ 0*t, z .+ 0*t
    # Buffers for displacements:
    ux, uy, uz = similar(xt), similar(yt), similar(zt) 

    # Composable motions: they need to be run sequentially
    for motion in Iterators.filter(is_composable, motion.types)
        displacement_x!(ux, motion, xt, yt, zt, t)
        displacement_y!(uy, motion, xt, yt, zt, t)
        displacement_z!(uz, motion, xt, yt, zt, t)
        xt .+= ux
        yt .+= uy
        zt .+= uz
    end
    # Additive motions: these motions can be run in parallel
    for motion in Iterators.filter(!is_composable, motion.types)
        displacement_x!(ux, motion, x, y, z, t)
        displacement_y!(uy, motion, x, y, z, t)
        displacement_z!(uz, motion, x, y, z, t)
        xt .+= ux
        yt .+= uy
        zt .+= uz
    end
    return xt, yt, zt
end

function times(motion::SimpleMotion)
    nodes = reduce(vcat, [times(type) for type in motion.types])
    nodes = unique(sort(nodes))
    return nodes
end

# utils: Sort Motion 
include("simplemotion/_utils.jl")

# Simple Motion Types:
# Non-periodic types: defined by an initial time (t_start), an end time (t_end) and a displacement      
include("simplemotion/Translation.jl")
include("simplemotion/Rotation.jl")
include("simplemotion/HeartBeat.jl")
# Periodic types: defined by the period, the temporal symmetry and a displacement (amplitude)
include("simplemotion/PeriodicTranslation.jl")
include("simplemotion/PeriodicRotation.jl")
include("simplemotion/PeriodicHeartBeat.jl")


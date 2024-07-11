"""
    motion = SimpleMotion(types)

SimpleMotion model. It allows for the definition of motion by means of simple parameters.
The `SimpleMotion` struct is composed by only one field, called `types`, 
which is a tuple of simple motion types. This tuple will contain as many elements
as simple motions we want to combine.

# Arguments
- `types`: (`::Tuple{Vararg{<:SimpleMotionType{T}}}`) tuple of simple motion types

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
abstract type SimpleMotion{T<:Real} <: Motion{T} end

Base.getindex(motion::SimpleMotion, p::Union{AbstractRange, AbstractVector, Colon, Integer}) = motion
Base.view(motion::SimpleMotion, p::Union{AbstractRange, AbstractVector, Colon, Integer}) = motion

include("simplemotion/Translation.jl")
include("simplemotion/Rotation.jl")
include("simplemotion/HeartBeat.jl")
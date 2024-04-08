# ------ SimpleMotionTypes
abstract type SimpleMotionType{T<:Real} end

"""
Simple Motion

x = x + ux
y = y + uy
z = z + uz

"""
# -------- SimpleMotion
mutable struct SimpleMotion{T<:Real} <: MotionModel{T}
    types::Vector{<:SimpleMotionType{T}}
end

Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                    q::Union{AbstractRange,AbstractVector,Colon}) = motion

function get_spin_coords(motion::SimpleMotion{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real}
    xi, yi, zi = x, y, z
    composable_motions = motion.types[is_composable.(motion.types)]
    sort!(composable_motions, by= m->m.ti)
    for motion in composable_motions
        xi, yi, zi = position_x(motion, xi, yi, zi, t), position_y(motion, xi, yi, zi, t), position_z(motion, xi, yi, zi, t)
    end
    xt = xi .+ reduce(.+, map((type) -> !is_composable(type) ? displacement_x(type, x, y, z, t) : 0, motion.types))
    yt = yi .+ reduce(.+, map((type) -> !is_composable(type) ? displacement_y(type, x, y, z, t) : 0, motion.types))
    zt = zi .+ reduce(.+, map((type) -> !is_composable(type) ? displacement_z(type, x, y, z, t) : 0, motion.types))
    return xt, yt, zt
end

function get_times(motion::SimpleMotion)
    times = reduce(vcat, [get_time_nodes(type) for type in motion.types])
    times = unique(sort(times))
    return times
end

# --------- Simple Motion Types: -------------
# Non-periodic types: defined by an initial time (ti), a final time (tf) and a displacement      
include("simplemotion/Translation.jl")
include("simplemotion/Rotation.jl")
include("simplemotion/HeartBeat.jl")
                                    
# Periodic types: defined by the period, the temporal symmetry and a displacement (amplitude)
include("simplemotion/PeriodicTranslation.jl")
include("simplemotion/PeriodicRotation.jl")
include("simplemotion/PeriodicHeartBeat.jl")

normalize_time_triangular(t, period, asymmetry) = begin
    t_rise = period * asymmetry
    t_fall = period * (1 - asymmetry)
    t_relative = mod.(t, period)
    t_unit = ifelse.(t_relative .< t_rise, t_relative ./ t_rise, 1 .- (t_relative .- t_rise) ./ t_fall)
    return t_unit
end


"""
Idea for motion pipeline (simple motion composition)

positions = x, y, z
for type in motion.types
positions = get_coords(type, positions...)
end
return positions
"""
# +(x::SimpleMotionType, y::SimpleMotionType)


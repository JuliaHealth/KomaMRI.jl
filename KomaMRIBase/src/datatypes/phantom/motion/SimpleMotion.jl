# ------ SimpleMotionType
abstract type SimpleMotionType{T<:Real} end

"""
Simple Motion

x = x + ux
y = y + uy
z = z + uz

"""
# -------- SimpleMotion
mutable struct SimpleMotion{T<:Real} <: MotionModel
    types::AbstractVector{SimpleMotionType{T}}
end

Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                    q::Union{AbstractRange,AbstractVector,Colon}) = motion

function get_spin_coords(motion::SimpleMotion{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real}
    xt = x .+ reduce(.+, map((type) -> displacement_x(type, x, y, z, t), motion.types))
    yt = y .+ reduce(.+, map((type) -> displacement_y(type, x, y, z, t), motion.types))
    zt = z .+ reduce(.+, map((type) -> displacement_z(type, x, y, z, t), motion.types))
    return xt, yt, zt
end

function get_range(motion::SimpleMotion)
    mini = minimum([get_range(type)[1] for type in motion.types])
    maxi = maximum([get_range(type)[2] for type in motion.types])
    return mini, maxi
end

# --------- Simple Motion Types: -------------
# Non-periodic types: defined by an initial time (ti), a final time (tf) and a displacement      
include("simplemotion/Translation.jl")
include("simplemotion/Rotation.jl")
#TODO: strain: how do I define it?
                                    
# Periodic types: defined by the period, the temporal symmetry and a displacement (amplitude)





"""
Idea for motion pipeline (simple motion composition)

positions = x, y, z
for type in motion.types
positions = get_coords(type, positions...)
end
return positions
"""
# +(x::SimpleMotionType, y::SimpleMotionType)


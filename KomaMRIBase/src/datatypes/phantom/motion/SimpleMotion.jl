# ------ SimpleMotionType
abstract type SimpleMotionType{T<:Real} end

# -------- SimpleMotion
mutable struct SimpleMotion{S <: SimpleMotionType} <: MotionModel
    types::AbstractVector{S}
end

Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::SimpleMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                    q::Union{AbstractRange,AbstractVector,Colon}) = motion


function displacement_x(motion::SimpleMotion{S}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real, S<:SimpleMotionType{T}}
    return reduce(.+, map((type) -> displacement_x(type, x, y, z, t), motion.types))
end

function displacement_y(motion::SimpleMotion{S}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real, S<:SimpleMotionType{T}}
    return reduce(.+, map((type) -> displacement_y(type, x, y, z, t), motion.types))
end

function displacement_z(motion::SimpleMotion{S}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real, S<:SimpleMotionType{T}}
    return reduce(.+, map((type) -> displacement_z(type, x, y, z, t), motion.types))
end

# --------- Simple Motion Types:
include("simplemotion/Translation.jl")
# include("simplemotion/Rotation.jl")


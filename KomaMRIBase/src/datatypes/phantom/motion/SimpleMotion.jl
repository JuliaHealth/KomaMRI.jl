# ------ SimpleMotionType
abstract type SimpleMotionType{T<:Real} end

# -------- SimpleMotion
mutable struct SimpleMotion{S <: SimpleMotionType} <: MotionModel
    types::AbstractVector{S}
end

Base.getindex(motion::SimpleMotion, p::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractRange, q::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractVector) = motion


function ux(motion::SimpleMotion{S}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real, S<:SimpleMotionType{T}}
    return reduce(.+, map((type) -> ux(type, x, y, z, t), motion.types))
end

function uy(motion::SimpleMotion{S}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real, S<:SimpleMotionType{T}}
    return reduce(.+, map((type) -> uy(type, x, y, z, t), motion.types))
end

function uz(motion::SimpleMotion{S}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real, S<:SimpleMotionType{T}}
    return reduce(.+, map((type) -> uz(type, x, y, z, t), motion.types))
end

# --------- Simple Motion Types:
include("simplemotion/Translation.jl")
# include("simplemotion/Rotation.jl")


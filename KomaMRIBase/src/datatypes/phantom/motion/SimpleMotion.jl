# ------ SimpleMotionType
abstract type SimpleMotionType{T<:Real} end

# -------- SimpleMotion
mutable struct SimpleMotion{T<:Real} <: MotionModel
    type::SimpleMotionType{T}
    ux::Function
	uy::Function
	uz::Function
end

# IDEA:
# struct SimpleMotion{T <: SimpleMotionType}
#     type::T
# end

Base.getindex(motion::SimpleMotion, p::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractRange, q::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractVector) = motion


function get_positions(motion::SimpleMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real}
	xt = x.+ motion.ux(x, y, z, t)
    yt = y.+ motion.uy(x, y, z, t)
    zt = z.+ motion.uz(x, y, z, t)
      
    return xt, yt, zt
end


# --------- Simple Motion Types:
include("simplemotion/Translation.jl")
# include("simplemotion/Rotation.jl")




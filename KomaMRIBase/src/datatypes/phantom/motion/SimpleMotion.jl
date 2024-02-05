@with_kw mutable struct SimpleMotion <: MotionModel
    ux::Function = (x,y,z,t)->0
	uy::Function = (x,y,z,t)->0
	uz::Function = (x,y,z,t)->0
end


Base.getindex(motion::SimpleMotion, p::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractRange, q::AbstractRange) = motion
Base.getindex(motion::SimpleMotion, p::AbstractVector) = motion


function get_displacements(motion::SimpleMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}) where {T<:Real}
	xt = x.+ motion.ux(x, y, z, t')
    yt = y.+ motion.uy(x, y, z, t')
    zt = z.+ motion.uz(x, y, z, t')

    xt = size(xt, 2) == 1 ? vec(xt) : xt
    yt = size(yt, 2) == 1 ? vec(yt) : yt
    zt = size(zt, 2) == 1 ? vec(zt) : zt
        
    xt, yt, zt, nothing
end



# ------ SimpleMotionType

abstract type SimpleMotionType end

# Simple Motion Types:
include("simplemotion/Rotation.jl")

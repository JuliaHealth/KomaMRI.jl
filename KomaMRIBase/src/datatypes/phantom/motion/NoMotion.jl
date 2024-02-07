@with_kw mutable struct NoMotion{T<:Real} <: MotionModel
    Δx::T = 0.0
    Δy::T = 0.0
    Δz::T = 0.0
end

Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                q::Union{AbstractRange,AbstractVector,Colon}) = motion

function get_positions(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real}
	x, y, z, nothing
end
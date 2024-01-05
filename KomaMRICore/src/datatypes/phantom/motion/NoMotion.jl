@with_kw mutable struct NoMotion{T} <: MotionModel{T}
    Δx::T = 0.0
    Δy::T = 0.0
    Δz::T = 0.0
end

Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}) = motion
Base.getindex(motion::NoMotion, p::Union{AbstractRange,AbstractVector,Colon}, 
                                q::Union{AbstractRange,AbstractVector,Colon}) = motion

function get_displacements(motion::NoMotion, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractVector{T}) where {T<:Real}
	0, 0, 0, nothing
end
struct NoMotion{T<:Real} <: AbstractMotion{T} end

Base.getindex(mv::NoMotion, p::Union{AbstractRange, AbstractVector, Colon, Integer}) = mv 
Base.view(mv::NoMotion, p::Union{AbstractRange, AbstractVector, Colon, Integer})     = mv

""" Addition of NoMotions """
Base.vcat(m1::NoMotion{T}, m2::AbstractMotion{T}, Ns1::Int, Ns2::Int) where {T<:Real} = m2
Base.vcat(m1::AbstractMotion{T}, m2::NoMotion{T}, Ns1::Int, Ns2::Int) where {T<:Real} = m1

Base.:(==)(m1::NoMotion{T}, m2::NoMotion{T}) where {T<:Real} = true
Base.:(≈)(m1::NoMotion{T}, m2::NoMotion{T}) where {T<:Real}  = true

function get_spin_coords(
    mv::NoMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T}
) where {T<:Real}
    return x, y, z
end

"""
    times = times(motion)
"""
times(mv::NoMotion{T}) where {T<:Real} = [zero(T)]

"""
    sort_motions!
"""
sort_motions!(mv::NoMotion) =  nothing
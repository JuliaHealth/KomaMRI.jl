const Interpolator1D = Interpolations.GriddedInterpolation{
    T,1,V,Itp,K
} where {
    T<:Real,
    V<:AbstractArray{<:Union{T,Bool}},
    Itp<:Interpolations.Gridded,
    K<:Tuple{AbstractVector{T}},
}

const Interpolator2D = Interpolations.GriddedInterpolation{
    T,2,V,Itp,K
} where {
    T<:Real,
    V<:AbstractArray{<:Union{T,Bool}},
    Itp<:Interpolations.Gridded,
    K<:Tuple{AbstractVector{T}, AbstractVector{T}},
}

"""
    ArbitraryMotion
"""
abstract type ArbitraryMotion{T<:Real} <: AbstractMotion{T} end

function Base.getindex(motion::ArbitraryMotion, p::Union{AbstractRange, AbstractVector, Colon, Integer})
    return typeof(motion)(motion.time, [getfield(motion, d)[p,:] for d in filter(x -> x != :time, fieldnames(typeof(motion)))]...)
end
function Base.view(motion::ArbitraryMotion, p::Union{AbstractRange, AbstractVector, Colon, Integer})
    return typeof(motion)(motion.time, [@view(getfield(motion, d)[p,:]) for d in filter(x -> x != :time, fieldnames(typeof(motion)))]...)
end

Base.:(==)(m1::ArbitraryMotion, m2::ArbitraryMotion) = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(typeof(m1))])
Base.:(≈)(m1::ArbitraryMotion, m2::ArbitraryMotion)  = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field)  ≈ getfield(m2, field) for field in fieldnames(typeof(m1))])

function GriddedInterpolation(nodes, A, ITP)
    return Interpolations.GriddedInterpolation{eltype(A), length(nodes), typeof(A), typeof(ITP), typeof(nodes)}(nodes, A, ITP)
end

function interpolate(d::AbstractArray{T}, ITPType, Ns::Val{1}) where {T<:Real}
    _, Nt = size(d)
    t = similar(d, Nt); copyto!(t, collect(range(zero(T), oneunit(T), Nt)))
    return GriddedInterpolation((t, ), d[:], ITPType)
end

function interpolate(d::AbstractArray{T}, ITPType, Ns::Val) where {T<:Real}
    Ns, Nt = size(d)
    id = similar(d, Ns); copyto!(id, collect(range(oneunit(T), T(Ns), Ns)))
    t = similar(d, Nt);  copyto!(t, collect(range(zero(T), oneunit(T), Nt)))
    return GriddedInterpolation((id, t), d, ITPType)
end

function resample(itp::Interpolator1D{T}, t::AbstractArray{T}) where {T<:Real}
    return itp.(t)
end

function resample(itp::Interpolator2D{T}, t::AbstractArray{T}) where {T<:Real}
    Ns = size(itp.coefs, 1)
    id = similar(itp.coefs, Ns)
    copyto!(id, collect(range(oneunit(T), T(Ns), Ns)))
    return itp.(id, t)
end

function displacement_x!(
    ux::AbstractArray{T},
    motion::ArbitraryMotion{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(motion.dx, Gridded(Linear()), Val(size(x,1)))
    ux .= resample(itp, unit_time(t, motion.time))
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
    motion::ArbitraryMotion{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(motion.dy, Gridded(Linear()), Val(size(x,1)))
    uy .= resample(itp, unit_time(t, motion.time))
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    motion::ArbitraryMotion{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(motion.dz, Gridded(Linear()), Val(size(x,1)))
    uz .= resample(itp, unit_time(t, motion.time))
    return nothing
end

include("arbitrarymotion/Trajectory.jl")
include("arbitrarymotion/FlowTrajectory.jl")
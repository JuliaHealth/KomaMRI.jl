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

abstract type ArbitraryAction{T<:Real} <: AbstractActionSpan{T} end

function Base.getindex(action::ArbitraryAction, p::Union{AbstractVector, Colon})
    return typeof(action)([getfield(action, d)[p,:] for d in fieldnames(typeof(action))]...)
end
function Base.view(action::ArbitraryAction, p::Union{AbstractVector, Colon})
    return typeof(action)([@view(getfield(action, d)[p,:]) for d in fieldnames(typeof(action))]...)
end

Base.:(==)(m1::ArbitraryAction, m2::ArbitraryAction) = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(typeof(m1))])
Base.:(≈)(m1::ArbitraryAction, m2::ArbitraryAction)  = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field)  ≈ getfield(m2, field) for field in fieldnames(typeof(m1))])

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
    # println(@view(t[1:3]))
    r = itp.(id, t)
    # println(@view(r[1:3]), '\n')
    return r
end

function displacement_x!(
    ux::AbstractArray{T},
    action::ArbitraryAction{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(action.dx, Gridded(Linear()), Val(size(action.dx,1)))
    ux .= resample(itp, t)
    return nothing
end

function displacement_y!(
    uy::AbstractArray{T},
    action::ArbitraryAction{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(action.dy, Gridded(Linear()), Val(size(action.dy,1)))
    uy .= resample(itp, t)
    return nothing
end

function displacement_z!(
    uz::AbstractArray{T},
    action::ArbitraryAction{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(action.dz, Gridded(Linear()), Val(size(action.dz,1)))
    uz .= resample(itp, t)
    return nothing
end




function displacement_x(
    action::ArbitraryAction{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(action.dx, Gridded(Linear()), Val(size(action.dx,1)))
    return resample(itp, t)
end

function displacement_y(
    action::ArbitraryAction{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(action.dy, Gridded(Linear()), Val(size(action.dy,1)))
    uy = resample(itp, t)
    m = minimum([size(uy,2), 8])
    # println("t:  ", @view(t[1, 1:m]))
    println("uy: ", @view(uy[1, 1:m]))
    print("\n")
    return uy
end

function displacement_z(
    action::ArbitraryAction{T},
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T},
    t::AbstractArray{T},
) where {T<:Real}
    itp = interpolate(action.dz, Gridded(Linear()), Val(size(action.dz,1)))
    return resample(itp, t)
end

include("arbitraryactions/Path.jl")
include("arbitraryactions/FlowPath.jl")
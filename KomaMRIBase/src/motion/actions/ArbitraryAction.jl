# We defined two types of Interpolation objects: Interpolator1D and Interpolator2D
# 1D is for interpolating for 1 spin
# 2D is for interpolating for 2 or more spins
# This dispatch based on the number of spins wouldn't be necessary if it weren't for this:
# https://github.com/JuliaMath/Interpolations.jl/issues/603
# 
# Once this issue is solved, this file should be simpler. 
# We should then be able to define a single method for functions:
#   - interpolate
#   - resample
# and delete the Interpolator1D and Interpolator2D definitions

const Interpolator1D = Interpolations.GriddedInterpolation{
    TCoefs,1,V,Itp,K
} where {
    TCoefs<:Real,
    TNodes<:Real,
    V<:AbstractArray{TCoefs},
    Itp<:Interpolations.Gridded,
    K<:Tuple{AbstractVector{TNodes}},
}

const Interpolator2D = Interpolations.GriddedInterpolation{
    TCoefs,2,V,Itp,K
} where {
    TCoefs<:Real,
    TNodes<:Real,
    V<:AbstractArray{TCoefs},
    Itp<:Interpolations.Gridded,
    K<:Tuple{AbstractVector{TNodes}, AbstractVector{TNodes}},
}

abstract type ArbitraryAction{T<:Real} <: AbstractAction{T} end

function Base.getindex(action::ArbitraryAction, p)
    return typeof(action)([getfield(action, d)[p,:] for d in fieldnames(typeof(action))]...)
end
function Base.view(action::ArbitraryAction, p)
    return typeof(action)([@view(getfield(action, d)[p,:]) for d in fieldnames(typeof(action))]...)
end

Base.:(==)(m1::ArbitraryAction, m2::ArbitraryAction) = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(typeof(m1))])
Base.:(≈)(m1::ArbitraryAction, m2::ArbitraryAction)  = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field)  ≈ getfield(m2, field) for field in fieldnames(typeof(m1))])

function GriddedInterpolation(nodes, A, ITP)
    return Interpolations.GriddedInterpolation{eltype(A), length(nodes), typeof(A), typeof(ITP), typeof(nodes)}(nodes, A, ITP)
end

function interpolate(d, ITPType, Ns::Val{1}, t)
    _, Nt = size(d)
    t_knots = _similar(t, Nt); copyto!(t_knots, collect(range(zero(eltype(t)), oneunit(eltype(t)), Nt)))
    return GriddedInterpolation((t_knots, ), d[:], ITPType)
end

function interpolate(d, ITPType, Ns::Val, t)
    Ns, Nt = size(d)
    id_knots = _similar(t, Ns); copyto!(id_knots, collect(range(oneunit(eltype(t)), eltype(t)(Ns), Ns)))
    t_knots  = _similar(t, Nt); copyto!(t_knots,  collect(range(zero(eltype(t)), oneunit(eltype(t)), Nt)))
    return GriddedInterpolation((id_knots, t_knots), d, ITPType)
end

function resample(itp::Interpolator1D, t)
    return itp.(t)
end

function resample(itp::Interpolator2D, t)
    Ns = size(itp.coefs, 1)
    id = _similar(t, Ns)
    copyto!(id, collect(range(oneunit(eltype(t)), eltype(t)(Ns), Ns)))
    return itp.(id, t)
end

function displacement_x!(ux, action::ArbitraryAction, x, y, z, t)
    itp = interpolate(action.dx, Gridded(Linear()), Val(size(action.dx,1)), t)
    ux .= resample(itp, t)
    return nothing
end

function displacement_y!(uy, action::ArbitraryAction, x, y, z, t)
    itp = interpolate(action.dy, Gridded(Linear()), Val(size(action.dy,1)), t)
    uy .= resample(itp, t)
    return nothing
end

function displacement_z!(uz, action::ArbitraryAction, x, y, z, t)
    itp = interpolate(action.dz, Gridded(Linear()), Val(size(action.dz,1)), t)
    uz .= resample(itp, t)
    return nothing
end

_similar(a, N) = similar(a, N)
_similar(a::Real, N) = zeros(typeof(a), N)

include("arbitraryactions/Path.jl")
include("arbitraryactions/FlowPath.jl")
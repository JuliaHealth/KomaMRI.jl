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

function interpolate_times(t, t_unit, periodic, tq)
    itp = GriddedInterpolation((t, ), t_unit, Gridded(Linear()))
    return extrapolate(itp, periodic ? Interpolations.Periodic() : Flat()).(tq)
end

_similar(a, N) = similar(a, N)
_similar(a::Real, N) = zeros(typeof(a), N)
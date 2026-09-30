# Keep the interpolate/resample interface used by older KomaMRICore versions.

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
    Itp<:Tuple{Interpolations.NoInterp, Interpolations.Gridded},
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
    return GriddedInterpolation((id_knots, t_knots), d, (Interpolations.NoInterp(), ITPType))
end

function resample(itp::Interpolator1D, t)
    return itp.(t)
end

function resample(itp::Interpolator2D, t)
    return itp.(itp.knots[1], t)
end

trajectory_knots(Nt, ::Type{T}) where {T<:Real} = range(zero(T), oneunit(T); length=Nt)

function trajectory_knots(Nt, t::AbstractArray)
    knots = similar(t, Nt)
    copyto!(knots, collect(trajectory_knots(Nt, eltype(t))))
    return knots
end

@inline function trajectory_point(knots, t)
    k = clamp(searchsortedfirst(knots, t) - 1, 1, length(knots) - 1)
    @inbounds w = (t - knots[k]) / (knots[k + 1] - knots[k])
    return k, w
end

@inline function trajectory_sample(d, spin, (k, w))
    @inbounds return (oneunit(w) - w) * d[spin, k] + w * d[spin, k + 1]
end

function resample_linear!(out, d, t::Real)
    point = trajectory_point(trajectory_knots(size(d, 2), typeof(t)), t)
    out .= trajectory_sample.(Ref(d), axes(d, 1), Ref(point))
    return nothing
end

function resample_linear!(out, d, t::AbstractArray)
    knots = trajectory_knots(size(d, 2), t)
    points = reshape(trajectory_point.(Ref(knots), t), 1, :)
    out .= trajectory_sample.(Ref(d), axes(d, 1), points)
    return nothing
end

@inline trajectory_next_column(knots, t) = min(max(searchsortedfirst(knots, t) - 1, 1) + 1, length(knots))

next_columns(d, t::Real) = @view d[:, trajectory_next_column(trajectory_knots(size(d, 2), typeof(t)), t)]

function next_columns(d, t::AbstractArray)
    knots = trajectory_knots(size(d, 2), t)
    return view(d, :, vec(trajectory_next_column.(Ref(knots), t)))
end

function interpolate_times(t, t_unit, periodic, tq)
    itp = GriddedInterpolation((t, ), t_unit, Gridded(Linear()))
    return extrapolate(itp, periodic ? Interpolations.Periodic() : Flat()).(tq)
end

_similar(a, N) = similar(a, N)
_similar(a::Real, N) = zeros(typeof(a), N)
# TODO: Consider different Extrapolations apart from periodic LinerInterpolator{T,ETPType}
#       Interpolator{T,Degree,ETPType}, 
#           Degree = Linear,Cubic.... 
#           ETPType = Periodic, Flat...

const Interpolator = Interpolations.GriddedInterpolation{
    T,N,V,Itp,K
} where {
    T<:Real,
    N,
    V<:AbstractArray{T},
    Itp<:Tuple{ Vararg{ Union{Interpolations.Gridded{Linear{Throw{OnGrid}}},Interpolations.NoInterp} } },
    K<:Tuple{Vararg{AbstractVector{T}}},
}

function GriddedInterpolation(
    nodes::NType,
    A::AType
) where {T<:Real, AType<:AbstractArray{T}, NType<:Tuple{Vararg{AbstractVector{T}}}}
    Ns, _ = size(A)
    if Ns > 1
        ITPType = Tuple{NoInterp, Gridded{Linear{Throw{OnGrid}}}}
        return Interpolations.GriddedInterpolation{T, 2, typeof(A), ITPType, typeof(nodes)}(nodes, A, (NoInterp(), Gridded(Linear())))
    else
        ITPType = Tuple{Gridded{Linear{Throw{OnGrid}}}}
        return Interpolations.GriddedInterpolation{T, 1, typeof(A[:]), ITPType, typeof((nodes[2], ))}((nodes[2], ), A[:], (Gridded(Linear()), ))
    end
end

"""
    motion = ArbitraryMotion(period_durations, dx, dy, dz)

ArbitraryMotion model. For this motion model, it is necessary to define 
motion for each spin independently, in x (`dx`), y (`dy`) and z (`dz`).
`dx`, `dy` and `dz` are three matrixes, of (``N_{spins}`` x ``N_{discrete\\,times}``) each.
This means that each row corresponds to a spin trajectory over a set of discrete time instants.
`period_durations` is a vector that contains the period for periodic (one element) or 
pseudo-periodic (two or more elements) motion.
The discrete time instants are calculated diving `period_durations` by ``N_{discrete\\,times}``.

This motion model is useful for defining arbitrarly complex motion, specially
for importing the spin trajectories from another source, like XCAT or a CFD.

# Arguments
- `period_durations`: (`Vector{T}`) 
- `dx`: (`::Array{T,2}`) matrix for displacements in x
- `dy`: (`::Array{T,2}`) matrix for displacements in y
- `dz`: (`::Array{T,2}`) matrix for displacements in z

# Returns
- `motion`: (`::ArbitraryMotion`) ArbitraryMotion struct

# Examples
```julia-repl
julia> motion = ArbitraryMotion(
            [1.0], 
            0.01.*rand(1000, 10), 
            0.01.*rand(1000, 10), 
            0.01.*rand(1000, 10)
        )
```
"""
struct ArbitraryMotion{T} <: MotionModel{T}
    t_start::T
    t_end::T
    dx::AbstractArray{T}
    dy::AbstractArray{T}
    dz::AbstractArray{T}
end

function Base.getindex(
    motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}
)
    return ArbitraryMotion(motion.t_start, motion.t_end, motion.dx[p,:], motion.dy[p,:], motion.dz[p,:])
end

Base.:(==)(m1::ArbitraryMotion, m2::ArbitraryMotion) = reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(ArbitraryMotion)])
Base.:(≈)(m1::ArbitraryMotion, m2::ArbitraryMotion)  = reduce(&, [getfield(m1, field) ≈ getfield(m2, field) for field in fieldnames(ArbitraryMotion)])

function Base.vcat(m1::ArbitraryMotion, m2::ArbitraryMotion)
    fields = []
    @assert (m1.t_start == m2.t_start) && (m1.t_end == m2.t_end) "starting and ending times must be the same"
    for field in (:dx, :dy, :dz)
        push!(fields, [getfield(m1, field); getfield(m2, field)])
    end
    return ArbitraryMotion(m1.t_start, m1.t_end, fields...)
end

"""
    limits = times(obj.motion)
"""
function times(motion::ArbitraryMotion)
    return range(motion.t_start, motion.t_end, length=size(motion.dx, 2))
end


function get_itp_functions(
    dx::AbstractArray{T},
    dy::AbstractArray{T},
    dz::AbstractArray{T}
    ) where {T<:Real}
    Ns, Nt = size(dx)

    t = similar(dx, Nt)
    t .= range(0,1,Nt)

    id = similar(dx, Ns)
    id .= 1:Ns

    itpx = GriddedInterpolation((id, t), dx)
    itpy = GriddedInterpolation((id, t), dy)
    itpz = GriddedInterpolation((id, t), dz)

    return itpx, itpy, itpz
end

function get_itp_results(
    itpx::Interpolator{T}, 
    itpy::Interpolator{T}, 
    itpz::Interpolator{T}, 
    t::AbstractArray{T}
) where {T<:Real}
    Ns = ndims(itpx.coefs) == 1 ? 1 : size(itpx.coefs,1)
    if Ns > 1
        # Grid
        idx = 1*id .+ 0*t # spin id
        t   = 0*id .+ 1*t # time instants
        return itpx.(idx, t), itpy.(idx, t), itpz.(idx, t)
    else
        return itpx.(t), itpy.(t), itpz.(t)
    end
end

function get_spin_coords(
    motion::ArbitraryMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T}
) where {T<:Real}
    itp = get_itp_functions(motion.dx, motion.dy, motion.dz)
    ux, uy, uz = get_itp_results(itp..., unit_time(t, motion.t_start, motion.t_end))
    return x .+ ux, y .+ uy, z .+ uz
end
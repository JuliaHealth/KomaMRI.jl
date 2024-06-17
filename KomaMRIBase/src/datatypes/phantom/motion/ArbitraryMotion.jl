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
    K<:Tuple{Vararg{AbstractRange{T}}},
}

function GriddedInterpolation(
    nodes::NType,
    A::AType
) where {T<:Real, AType<:AbstractArray{T}, NType<:Tuple{Vararg{AbstractRange{T}}}}
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
    @assert (m1.t_start == m2.t_start) && (m1.t_end == m2.t_end) "Starting and ending times must be the same"
    return ArbitraryMotion(m1.t_start, m1.t_end, [m1.dx; m2.dx], [m1.dy; m2.dy], [m1.dz; m2.dz])
end

"""
    limits = times(obj.motion)
"""
function times(motion::ArbitraryMotion)
    return range(motion.t_start, motion.t_end, length=size(motion.dx, 2))
end


function interpolate(motion::ArbitraryMotion{T}) where {T<:Real}
    Ns, Nt = size(motion.dx)
    itpx = GriddedInterpolation((one(T):Ns, range(0,one(T),Nt)), motion.dx)
    itpy = GriddedInterpolation((one(T):Ns, range(0,one(T),Nt)), motion.dy)
    itpz = GriddedInterpolation((one(T):Ns, range(0,one(T),Nt)), motion.dz)
    return itpx, itpy, itpz
end

function resample(
    itpx::Interpolator{T}, 
    itpy::Interpolator{T}, 
    itpz::Interpolator{T}, 
    t::AbstractArray{T}
) where {T<:Real}
    Ns = ndims(itpx.coefs) == 1 ? 1 : size(itpx.coefs,1)
    if Ns > 1
        return itpx.(1:Ns, t), itpy.(1:Ns, t), itpz.(1:Ns, t)
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
    motion_functions = interpolate(motion)
    ux, uy, uz = resample(motion_functions..., unit_time(t, motion.t_start, motion.t_end))
    return x .+ ux, y .+ uy, z .+ uz
end
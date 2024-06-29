# TODO: Consider different Extrapolations apart from periodic LinerInterpolator{T,ETPType}
#       Interpolator{T,Degree,ETPType}, 
#           Degree = Linear,Cubic.... 
#           ETPType = Periodic, Flat...

const Interpolator1D = Interpolations.GriddedInterpolation{
    T,1,V,Itp,K
} where {
    T<:Real,
    V<:AbstractArray{T},
    Itp<:Interpolations.Gridded{Linear{Throw{OnGrid}}},
    K<:Tuple{AbstractVector{T}},
}

const Interpolator2D = Interpolations.GriddedInterpolation{
    T,2,V,Itp,K
} where {
    T<:Real,
    V<:AbstractArray{T},
    Itp<:Interpolations.Gridded{Linear{Throw{OnGrid}}},
    K<:Tuple{AbstractVector{T}, AbstractVector{T}},
}

"""
    motion = ArbitraryMotion(t_start, t_end, dx, dy, dz)

ArbitraryMotion model. For this motion model, it is necessary to define 
motion for each spin independently, in x (`dx`), y (`dy`) and z (`dz`).
`dx`, `dy` and `dz` are three matrixes, of (``N_{spins}`` x ``N_{discrete\\,times}``) each.
This means that each row corresponds to a spin trajectory over a set of discrete time instants.
The motion duration is determined by `t_start` and `t_end`.
The discrete time instants are evenly spaced, (``dt = \frac{t_{end} - t_{start}}{N_{discrete\\,times}}``). 

This motion model is useful for defining arbitrarly complex motion, specially
for importing the spin trajectories from another source, like XCAT or a CFD.

# Arguments
- `t_start`: (`::T`, `[s]`) 
- `t_end`: (`::T`, `[s]`) 
- `dx`: (`::Array{T,2}`, `[m]`) matrix for displacements in x
- `dy`: (`::Array{T,2}`, `[m]`) matrix for displacements in y
- `dz`: (`::Array{T,2}`, `[m]`) matrix for displacements in z

# Returns
- `motion`: (`::ArbitraryMotion`) ArbitraryMotion struct

# Examples
```julia-repl
julia> motion = ArbitraryMotion(
            0.0,
            1.0, 
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
function Base.view(
    motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}
)
    return ArbitraryMotion(motion.t_start, motion.t_end, @view(motion.dx[p,:]), @view(motion.dy[p,:]), @view(motion.dz[p,:]))
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

function GriddedInterpolation(nodes, A, ITP)
    return Interpolations.GriddedInterpolation{eltype(A), length(nodes), typeof(A), typeof(ITP), typeof(nodes)}(nodes, A, ITP)
end

function interpolate(motion::ArbitraryMotion{T}, Ns::Val{1}) where {T<:Real}
    _, Nt = size(motion.dx)
    t = similar(motion.dx, Nt); copyto!(t, collect(range(zero(T), oneunit(T), Nt)))
    itpx = GriddedInterpolation((t, ), motion.dx[:], Gridded(Linear()))
    itpy = GriddedInterpolation((t, ), motion.dy[:], Gridded(Linear()))
    itpz = GriddedInterpolation((t, ), motion.dz[:], Gridded(Linear()))
    return itpx, itpy, itpz
end

function interpolate(motion::ArbitraryMotion{T}, Ns::Val) where {T<:Real}
    Ns, Nt = size(motion.dx)
    id = similar(motion.dx, Ns); copyto!(id, collect(range(oneunit(T), T(Ns), Ns)))
    t = similar(motion.dx, Nt);  copyto!(t, collect(range(zero(T), oneunit(T), Nt)))
    itpx = GriddedInterpolation((id, t), motion.dx, Gridded(Linear()))
    itpy = GriddedInterpolation((id, t), motion.dy, Gridded(Linear()))
    itpz = GriddedInterpolation((id, t), motion.dz, Gridded(Linear()))
    return itpx, itpy, itpz
end

function resample(itpx::Interpolator1D{T}, itpy::Interpolator1D{T}, itpz::Interpolator1D{T}, t::AbstractArray{T}) where {T<:Real}
    return itpx.(t), itpy.(t), itpz.(t)
end

function resample(itpx::Interpolator2D{T}, itpy::Interpolator2D{T}, itpz::Interpolator2D{T}, t::AbstractArray{T}) where {T<:Real}
    Ns = size(itpx.coefs, 1)
    id = similar(itpx.coefs, Ns)
    copyto!(id, collect(range(oneunit(T), T(Ns), Ns)))
    return itpx.(id, t), itpy.(id, t), itpz.(id, t)
end

function get_spin_coords(
    motion::ArbitraryMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T}
) where {T<:Real}
    motion_functions = interpolate(motion, Val(size(x,1)))
    ux, uy, uz = resample(motion_functions..., unit_time(t, motion.t_start, motion.t_end))
    return x .+ ux, y .+ uy, z .+ uz
end
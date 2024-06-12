# TODO: Consider different Extrapolations apart from periodic LinerInterpolator{T,ETPType}
#       Interpolator{T,Degree,ETPType}, 
#           Degree = Linear,Cubic.... 
#           ETPType = Periodic, Flat...

const Interpolator = Interpolations.Extrapolation{
    T,
    N,
    Interpolations.GriddedInterpolation{
        T,
        N,
        V,
        Itp,
        K
    },
    Itp,
    Interpolations.Periodic{Nothing}
} where {
    T<:Real,
    N,
    V<:AbstractArray{T},
    K<:Tuple{Vararg{AbstractVector{T}}},
    Itp<:Tuple{Vararg{Union{Interpolations.Gridded{Linear{Throw{OnGrid}}}, Interpolations.NoInterp}}}
}

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
    period_durations::Vector{T}
    dx::Matrix{T}
    dy::Matrix{T}
    dz::Matrix{T}
end

function Base.getindex(
    motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}
)
    fields = []
    for field in fieldnames(ArbitraryMotion)
        if field in (:dx, :dy, :dz)
            push!(fields, getfield(motion, field)[p, :])
        else
            push!(fields, getfield(motion, field))
        end
    end
    return ArbitraryMotion(fields...)
end

Base.:(==)(m1::ArbitraryMotion, m2::ArbitraryMotion) = reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(ArbitraryMotion)])
Base.:(≈)(m1::ArbitraryMotion, m2::ArbitraryMotion)  = reduce(&, [getfield(m1, field) ≈ getfield(m2, field) for field in fieldnames(ArbitraryMotion)])

function Base.vcat(m1::ArbitraryMotion, m2::ArbitraryMotion)
    fields = []
    @assert m1.period_durations == m2.period_durations "period_durations of both ArbitraryMotions must be the same"
    for field in
        Iterators.filter(x -> !(x == :period_durations), fieldnames(ArbitraryMotion))
        push!(fields, [getfield(m1, field); getfield(m2, field)])
    end
    return ArbitraryMotion(m1.period_durations, fields...)
end

"""
    limits = times(obj.motion)
"""
function times(motion::ArbitraryMotion)
    period_durations = motion.period_durations
    num_pieces = size(motion.dx)[2] + 1
    return times(period_durations, num_pieces)
end

function times(period_durations::AbstractVector, num_pieces::Int)
    # Pre-allocating memory
    limits = zeros(eltype(period_durations), num_pieces * length(period_durations) + 1)

    idx = 1
    for i in 1:length(period_durations)
        segment_increment = period_durations[i] / num_pieces
        cumulative_sum = limits[idx]  # Start from the last computed value in limits
        for j in 1:num_pieces
            cumulative_sum += segment_increment
            limits[idx + 1] = cumulative_sum
            idx += 1
        end
    end
    return limits
end


function get_itp_functions(motion::ArbitraryMotion{T}) where {T<:Real}
    Ns = size(motion.dx)[1]
    dx = hcat(repeat(hcat(zeros(T ,Ns, 1), motion.dx), 1, length(motion.period_durations)), zeros(T ,Ns, 1))
    dy = hcat(repeat(hcat(zeros(T ,Ns, 1), motion.dy), 1, length(motion.period_durations)), zeros(T ,Ns, 1))
    dz = hcat(repeat(hcat(zeros(T ,Ns, 1), motion.dz), 1, length(motion.period_durations)), zeros(T ,Ns, 1))
    if Ns > 1
        nodes = ([i*one(T) for i=1:Ns], times(motion))
        itpx = extrapolate(interpolate(nodes, dx, (NoInterp(), Gridded(Linear()))), Periodic())
        itpy = extrapolate(interpolate(nodes, dy, (NoInterp(), Gridded(Linear()))), Periodic())
        itpz = extrapolate(interpolate(nodes, dz, (NoInterp(), Gridded(Linear()))), Periodic())
    else
        nodes = (times(motion), )
        itpx = extrapolate(interpolate(nodes, dx[:], (Gridded(Linear()), )), Periodic())
        itpy = extrapolate(interpolate(nodes, dy[:], (Gridded(Linear()), )), Periodic())
        itpz = extrapolate(interpolate(nodes, dz[:], (Gridded(Linear()), )), Periodic())
    end
    return itpx, itpy, itpz
end

function get_itp_results(
    itpx::Interpolator{T}, 
    itpy::Interpolator{T}, 
    itpz::Interpolator{T}, 
    t::AbstractArray{T}, 
    range::AbstractRange
) where {T<:Real}
    Ns = length(range)
    if Ns > 1
        id = similar(t, Ns)
        id .= range
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
    x::Vector{T},
    y::Vector{T},
    z::Vector{T},
    t::AbstractArray{T}
) where {T<:Real}
    itp = get_itp_functions(motion)
    ux, uy, uz = get_itp_results(itp..., t, 1:length(x))
    return x .+ ux, y .+ uy, z .+ uz
end

function initialize_motion(motion::ArbitraryMotion)
    itp = get_itp_functions(motion)
    return ExplicitArbitraryMotion(itp..., 1:size(motion.dx)[1])
end


"""
    motion = ExplicitArbitraryMotion(period_durations, dx, dy, dz)

ExplicitArbitraryMotion model.

"""
mutable struct ExplicitArbitraryMotion{T} <: MotionModel{T}
    itpx::Interpolator{T}
    itpy::Interpolator{T}
    itpz::Interpolator{T}
    range::Union{AbstractRange,AbstractVector{T},Colon}
end

function Base.getindex(
    motion::ExplicitArbitraryMotion{T}, p::Union{AbstractRange,AbstractVector{T},Colon}
) where {T<:Real}
    motion.range = p
    return motion
end
  

function get_spin_coords(
    motion::ExplicitArbitraryMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T}
) where {T<:Real}
    ux, uy, uz = get_itp_results(motion.itpx, motion.itpy, motion.itpz, t, motion.range)
    return x .+ ux, y .+ uy, z .+ uz
end
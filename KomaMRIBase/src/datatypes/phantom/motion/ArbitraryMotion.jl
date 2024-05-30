# TODO: Consider different Extrapolations apart from periodic LinerInterpolator{T,ETPType}
#       Interpolator{T,Degree,ETPType}, 
#           Degree = Linear,Cubic.... 
#           ETPType = Periodic, Flat...
const LinearInterpolator = Interpolations.Extrapolation{
    T,
    1,
    Interpolations.GriddedInterpolation{T,1,V,Gridded{Linear{Throw{OnGrid}}},Tuple{V}},
    Gridded{Linear{Throw{OnGrid}}},
    Periodic{Nothing},
} where {T<:Real,V<:AbstractVector{T}}

"""
    motion = ArbitraryMotion(period_durations, dx, dy, dz)

ArbitraryMotion model. For this motion model, it is necessary to define 
motion for each spin independently, in x (`dx`), y (`dy`) and z (`dz`).
`dx`, `dy` and `dz` are three matrixes, of (``N_{spins}`` x ``N_{discrete\,times}``) each.
This means that each row corresponds to a spin trajectory over a set of discrete time instants.
`period_durations` is a vector that contains the period for periodic (one element) or 
pseudo-periodic (two or more elements) motion.
The discrete time instants are calculated diving `period_durations` by ``N_{discrete\,times}``.

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
struct ArbitraryMotion{T<:Real,V<:AbstractVector{T}} <: MotionModel{T}
    period_durations::Vector{T}
    dx::Array{T,2}
    dy::Array{T,2}
    dz::Array{T,2}
    ux::Vector{LinearInterpolator{T,V}}
    uy::Vector{LinearInterpolator{T,V}}
    uz::Vector{LinearInterpolator{T,V}}
end

function ArbitraryMotion(
    period_durations::AbstractVector{T},
    dx::AbstractArray{T,2},
    dy::AbstractArray{T,2},
    dz::AbstractArray{T,2},
) where {T<:Real}
    @warn "Note that ArbitraryMotion is under development so it is not optimized so far" maxlog = 1
    Ns = size(dx)[1]
    num_pieces = size(dx)[2] + 1
    limits = times(period_durations, num_pieces)

    #! format: off
    Δ = zeros(Ns,length(limits),4)
    Δ[:,:,1] = hcat(repeat(hcat(zeros(Ns,1),dx),1,length(period_durations)),zeros(Ns,1))
    Δ[:,:,2] = hcat(repeat(hcat(zeros(Ns,1),dy),1,length(period_durations)),zeros(Ns,1))
    Δ[:,:,3] = hcat(repeat(hcat(zeros(Ns,1),dz),1,length(period_durations)),zeros(Ns,1))
   
    etpx = [extrapolate(interpolate((limits,), Δ[i,:,1], Gridded(Linear())), Periodic()) for i in 1:Ns]
    etpy = [extrapolate(interpolate((limits,), Δ[i,:,2], Gridded(Linear())), Periodic()) for i in 1:Ns]
    etpz = [extrapolate(interpolate((limits,), Δ[i,:,3], Gridded(Linear())), Periodic()) for i in 1:Ns]
    #! format: on

    return ArbitraryMotion(period_durations, dx, dy, dz, etpx, etpy, etpz)
end

function Base.getindex(
    motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}
)
    fields = []
    for field in fieldnames(ArbitraryMotion)
        if field in (:dx, :dy, :dz)
            push!(fields, getfield(motion, field)[p, :])
        elseif field in (:ux, :uy, :uz)
            push!(fields, getfield(motion, field)[p])
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

# TODO: Calculate interpolation functions "on the fly"
function get_spin_coords(
    motion::ArbitraryMotion{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    z::AbstractVector{T},
    t::AbstractArray{T},
) where {T<:Real}
    xt = x .+ reduce(vcat, [etp.(t) for etp in motion.ux])
    yt = y .+ reduce(vcat, [etp.(t) for etp in motion.uy])
    zt = z .+ reduce(vcat, [etp.(t) for etp in motion.uz])
    return xt, yt, zt
end

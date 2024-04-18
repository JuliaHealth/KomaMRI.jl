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
Arbitrary Motion

x = x + ux
"""
struct ArbitraryMotion{T<:Real,V<:AbstractVector{T}} <: MotionModel{T}
    period_durations::AbstractVector{T}
    dx::AbstractArray{T,2}
    dy::AbstractArray{T,2}
    dz::AbstractArray{T,2}
    ux::Vector{LinearInterpolator{T,V}}
    uy::Vector{LinearInterpolator{T,V}}
    uz::Vector{LinearInterpolator{T,V}}
end

function ArbitraryMotion(
    period_durations::AbstractVector{T},
    Δx::AbstractArray{T,2},
    Δy::AbstractArray{T,2},
    Δz::AbstractArray{T,2},
) where {T<:Real}
    Ns = size(Δx)[1]
    num_pieces = size(Δx)[2] + 1
    limits = time_nodes(period_durations, num_pieces)

    #! format: off
    Δ = zeros(Ns,length(limits),4)
    Δ[:,:,1] = hcat(repeat(hcat(zeros(Ns,1),Δx),1,length(period_durations)),zeros(Ns,1))
    Δ[:,:,2] = hcat(repeat(hcat(zeros(Ns,1),Δy),1,length(period_durations)),zeros(Ns,1))
    Δ[:,:,3] = hcat(repeat(hcat(zeros(Ns,1),Δz),1,length(period_durations)),zeros(Ns,1))
   
    etpx = [extrapolate(interpolate((limits,), Δ[i,:,1], Gridded(Linear())), Periodic()) for i in 1:Ns]
    etpy = [extrapolate(interpolate((limits,), Δ[i,:,2], Gridded(Linear())), Periodic()) for i in 1:Ns]
    etpz = [extrapolate(interpolate((limits,), Δ[i,:,3], Gridded(Linear())), Periodic()) for i in 1:Ns]
    #! format: on

    return ArbitraryMotion(period_durations, Δx, Δy, Δz, etpx, etpy, etpz)
end

function Base.getindex(
    motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}
)
    fields = []
    for field in filter(x -> (x in [:dx, :dy, :dz]), [fieldnames(ArbitraryMotion)...])
        push!(fields, getfield(motion, field)[p, :])
    end
    for field in filter(x -> (x in [:ux, :uy, :uz]), [fieldnames(ArbitraryMotion)...])
        push!(fields, getfield(motion, field)[p])
    end
    return ArbitraryMotion(motion.period_durations, fields...)
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
    limits = time_nodes(obj.motion)
"""
# Revise this function to make it more efficient
function time_nodes(motion::ArbitraryMotion)
    period_durations = motion.period_durations
    num_pieces = size(motion.dx)[2] + 1
    return time_nodes(period_durations, num_pieces)
end

function time_nodes(period_durations::AbstractVector, num_pieces::Int)
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

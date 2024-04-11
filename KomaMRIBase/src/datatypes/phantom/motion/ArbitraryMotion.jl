# TODO: Consider different Extrapolations apart from periodic LinerInterpolator{T,ETPType}
#       Interpolator{T,Degree,ETPType}, 
#           Degree = Linear,Cubic.... 
#           ETPType = Periodic, Flat...
const LinearInterpolator = Interpolations.Extrapolation{T, 1, Interpolations.GriddedInterpolation{T, 1, V, Gridded{Linear{Throw{OnGrid}}}, Tuple{V}}, Gridded{Linear{Throw{OnGrid}}}, Periodic{Nothing}} where {T<:Real, V<:AbstractVector{T}}

"""
Arbitrary Motion

x = x + ux
"""
mutable struct ArbitraryMotion{T<:Real,V<:AbstractVector{T}} <: MotionModel{T}
    duration::AbstractVector{T}
    dx::AbstractArray{T,2}
    dy::AbstractArray{T,2}
    dz::AbstractArray{T,2}
    ux::Vector{LinearInterpolator{T,V}}
    uy::Vector{LinearInterpolator{T,V}}
    uz::Vector{LinearInterpolator{T,V}}
end

# TODO:
# mutable struct ArbitraryMotion{T<:Real} <: MotionModel
#     dur::AbstractVector{T}
#     Δx::AbstractArray{T, 2}
#     Δy::AbstractArray{T, 2}
#     Δz::AbstractArray{T, 2}
# end

# Optimize: ver https://github.com/cncastillo/KomaMRI.jl/issues/73
# t0 = [0; cumsum(dur)] 
# time = repeat(.., [1, length(dur)])
# time = time .+ t0'
# time = time[:]

function ArbitraryMotion( dur::AbstractVector{T},
                          Δx::AbstractArray{T, 2},
                          Δy::AbstractArray{T, 2},
                          Δz::AbstractArray{T, 2}) where {T<:Real}

    Ns = size(Δx)[1]
    K = size(Δx)[2] + 1
    limits = get_pieces_limits(dur,K)

    Δ = zeros(Ns,length(limits),4)
    Δ[:,:,1] = hcat(repeat(hcat(zeros(Ns,1),Δx),1,length(dur)),zeros(Ns,1))
    Δ[:,:,2] = hcat(repeat(hcat(zeros(Ns,1),Δy),1,length(dur)),zeros(Ns,1))
    Δ[:,:,3] = hcat(repeat(hcat(zeros(Ns,1),Δz),1,length(dur)),zeros(Ns,1))

    etpx = [extrapolate(interpolate((limits,), Δ[i,:,1], Gridded(Linear())), Periodic()) for i in 1:Ns]
    etpy = [extrapolate(interpolate((limits,), Δ[i,:,2], Gridded(Linear())), Periodic()) for i in 1:Ns]
    etpz = [extrapolate(interpolate((limits,), Δ[i,:,3], Gridded(Linear())), Periodic()) for i in 1:Ns]

    return ArbitraryMotion(dur, Δx, Δy, Δz, etpx, etpy, etpz)
end


"""
    limits = get_pieces_limits(obj.motion)

Returns the pieces limits from dur and K values

Example: -----------------------------
    motion.dur = [1, 0.5]
    motion.K = 4

    limits = [0, 0.25, 0.5, 0.75, 1, 1.125, 1.25, 1.375, 1.5]
--------------------------------------
"""
# Revise this function to make it more efficient
function get_pieces_limits(dur::AbstractVector, K::Int)
	steps = dur/K
	mat = reduce(hcat,[steps for i in 1:K])'
	limits = reshape(mat,(K*length(dur),))
	cumsum!(limits,limits)
	limits = vcat(0,limits)
    limits
end

Base.getindex(motion::ArbitraryMotion, p::Union{AbstractRange,AbstractVector,Colon}) = begin
    return ArbitraryMotion(
        motion.duration,
        motion.dx[p,:],
        motion.dy[p,:],
        motion.dz[p,:],
        motion.ux[p],
        motion.uy[p],
        motion.uz[p]
    )
end

Base.:(==)(m1::ArbitraryMotion, m2::ArbitraryMotion) = begin
    m1.duration == m2.duration  &&
    m1.dx       == m2.dx        &&
    m1.dy       == m2.dy        &&
    m1.dz       == m2.dz        &&
    m1.ux       == m2.ux        &&
    m1.uy       == m2.uy        &&
    m1.uz       == m2.uz
end

# TODO: Calculate interpolation functions "on the fly"
function get_spin_coords(motion::ArbitraryMotion{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t::AbstractArray{T}) where {T<:Real}
    xt = x .+ reduce(vcat, [etp.(t) for etp in motion.ux])
    yt = y .+ reduce(vcat, [etp.(t) for etp in motion.uy])
    zt = z .+ reduce(vcat, [etp.(t) for etp in motion.uz])
    return xt, yt, zt
end

function get_times(motion::ArbitraryMotion)
    return get_pieces_limits(motion.duration, size(motion.dx)[2] + 1)
end



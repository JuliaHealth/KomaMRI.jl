# TODO: Consider different Extrapolations apart from periodic LinerInterpolator{T,ETPType}
#       Interpolator{T,Degree,ETPType}, 
#           Degree = Linear,Cubic.... 
#           ETPType = Periodic, Flat...
"""
Arbitrary Motion

x = x + ux
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


function get_itp_functions(motion::ArbitraryMotion{T}, Ns::Int) where {T<:Real}
    dx = hcat(repeat(hcat(zeros(Ns, 1), motion.dx), 1, length(motion.period_durations)), zeros(Ns, 1))
    dy = hcat(repeat(hcat(zeros(Ns, 1), motion.dy), 1, length(motion.period_durations)), zeros(Ns, 1))
    dz = hcat(repeat(hcat(zeros(Ns, 1), motion.dz), 1, length(motion.period_durations)), zeros(Ns, 1))
    if Ns > 1
        nodes = (times(motion), [i for i=1:Ns])
        itpx = extrapolate(interpolate(nodes, dx, (Gridded(Linear(), NoInterp()))), Periodic())
        itpy = extrapolate(interpolate(nodes, dy, (Gridded(Linear(), NoInterp()))), Periodic())
        itpz = extrapolate(interpolate(nodes, dz, (Gridded(Linear(), NoInterp()))), Periodic())
    else
        nodes = (times(motion), )
        itpx = extrapolate(interpolate(nodes, dx[:], (Gridded(Linear()), )), Periodic())
        itpy = extrapolate(interpolate(nodes, dy[:], (Gridded(Linear()), )), Periodic())
        itpz = extrapolate(interpolate(nodes, dz[:], (Gridded(Linear()), )), Periodic())
    end
    return itpx, itpy, itpz
end

function get_itp_results(itpx, itpy, itpz, t, Ns)
    if Ns > 1
        id = 0 .* similar(t, Ns) .+ (1:Ns)
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
    t::AbstractArray{T},
) where {T<:Real}
    Ns = size(motion.dx)[1]
    itp = get_itp_functions(motion, Ns)
    ux, uy, uz = get_itp_results(itp..., t, Ns)
    return x .+ ux, y .+ uy, z .+ uz
end
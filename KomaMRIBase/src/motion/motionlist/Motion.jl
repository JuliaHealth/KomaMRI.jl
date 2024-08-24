"""
    m = Motion(action, time, spins)

Motion struct. (...)
"""
@with_kw mutable struct Motion{T<:Real}
    action::AbstractActionSpan{T}
    time::AbstractTimeSpan{T}
    spins::AbstractSpinSpan
end

""" Constructors """
function Motion(action::A, time::TS, spins::AbstractSpinSpan) where {T<:Real, A<:AbstractActionSpan{T}, TS<:AbstractTimeSpan{T}}
    return Motion{T}(action, time, spins)
end
function Motion(action::A, time::TS) where {T<:Real, A<:AbstractActionSpan{T}, TS<:AbstractTimeSpan{T}}
    return Motion{T}(action, time, AllSpins())
end
function Motion(action::A) where {T<:Real, A<:AbstractActionSpan{T}}
    return Motion{T}(action, TimeRange(zero(T)))
end
function Motion(action::A, time::TS, range::Colon) where {T<:Real, A<:AbstractActionSpan{T}, TS<:AbstractTimeSpan{T}}
    return Motion{T}(action, time, AllSpins())
end
function Motion(action::A, time::TS, range::AbstractVector) where {T<:Real, A<:AbstractActionSpan{T}, TS<:AbstractTimeSpan{T}}
    return Motion{T}(action, time, SpinRange(range))
end

# Custom constructors
function Translate(dx, dy, dz, time=TimeRange(0.0), spins=AllSpins())
    return Motion(Translate(dx, dy, dz), time, spins)
end
function Rotate(pitch, roll, yaw, time=TimeRange(0.0), spins=AllSpins())
    return Motion(Rotate(pitch, roll, yaw), time, spins)
end
function HeartBeat(circumferential_strain, radial_strain, longitudinal_strain, time=TimeRange(0.0), spins=AllSpins())
    return Motion(HeartBeat(circumferential_strain, radial_strain, longitudinal_strain), time, spins)
end
function Path(dx, dy, dz, time=TimeRange(0.0), spins=AllSpins())
    return Motion(Path(dx, dy, dz), time, spins)
end
function FlowPath(dx, dy, dz, spin_reset, time=TimeRange(0.0), spins=AllSpins())
    return Motion(FlowPath(dx, dy, dz, spin_reset), time, spins)
end

""" Compare two Motions """
Base.:(==)(m1::Motion, m2::Motion) = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(typeof(m1))])
Base.:(≈)(m1::Motion, m2::Motion)  = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field)  ≈ getfield(m2, field) for field in fieldnames(typeof(m1))])

""" Motion sub-group """
function Base.getindex(m::Motion, p::AbstractVector)
    idx, spin_range = m.spins[p]
    return Motion(m.action[idx], m.time, spin_range)
end
function Base.view(m::Motion, p::AbstractVector)
    idx, spin_range = @view(m.spins[p])
    return Motion(@view(m.action[idx]), m.time, spin_range)
end

# Auxiliary functions
times(m::Motion) = times(m.time)
add_motion!(motion_array, motion) = has_spins(motion.spins) ? push!(motion_array, motion) : nothing
is_composable(m::Motion) = is_composable(m.action)
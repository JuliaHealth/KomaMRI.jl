"""
    timecurve = TimeCurve(t, t_unit, periodic, duration)

TimeCurve struct. It is a specialized type that defines a time curve. 
(...)

# Arguments
- `t`: (`::AbstractVector{<:Real}`, `[s]`) time vector
- `t_unit`: (`::AbstractVector{<:Real}`) y vector, it needs to be scaled between 0 and 1
- `periodic`: (`::Bool`, `=false`) indicates whether the time curve should be periodically repeated
- `duration`: (`::Union{<:Real,AbstractVector{<:Real}}`, `=1.0`)

# Returns
- `timecurve`: (`::TimeCurve`) TimeCurve struct

# Examples
```julia-repl
julia> timecurve = TimeCurve(t=[0.0, 0.1, 0.3, 0.4], t_unit=[0.0, 0.6, 0.2, 0.0], periodic=true)
```
"""
@with_kw struct TimeCurve{T<:Real}
    t::AbstractVector{T}
    t_unit::AbstractVector{T}
    periodic::Bool                       = false
    duration::Union{T,AbstractVector{T}} = oneunit(eltype(t))
    t_start::T                           = t[1]
    t_end::T                             = t[end]
    @assert check_unique(t) "Vector t=$(t) contains duplicate elements. Please ensure all elements in t are unique and try again"
end

check_unique(t) = true
check_unique(t::Vector) = length(t) == length(unique(t))

# Main Constructors
TimeCurve(t, t_unit, periodic, duration) = TimeCurve(t=t, t_unit=t_unit, periodic=periodic, duration=duration)
TimeCurve(t, t_unit) = TimeCurve(t=t, t_unit=t_unit)
# Custom constructors
# --- TimeRange
TimeRange(t_start::T, t_end::T) where T = TimeCurve(t=[t_start, t_end], t_unit=[zero(T), oneunit(T)])
TimeRange(; t_start=0.0, t_end=1.0)     = TimeRange(t_start, t_end)
# --- Periodic
Periodic(period::T, asymmetry::T) where T = TimeCurve(t=[zero(T), period*asymmetry, period], t_unit=[zero(T), oneunit(T), zero(T)])
Periodic(; period=1.0, asymmetry=0.5)     = Periodic(period, asymmetry)

""" Compare two TimeCurves """
Base.:(==)(t1::TimeCurve, t2::TimeCurve) = reduce(&, [getfield(t1, field) == getfield(t2, field) for field in fieldnames(typeof(t1))])
Base.:(≈)(t1::TimeCurve, t2::TimeCurve)  = reduce(&, [getfield(t1, field)  ≈ getfield(t2, field) for field in fieldnames(typeof(t1))])

""" times """
function times(t, dur::AbstractVector)
    tr      = repeat(t, length(dur))
    scale   = repeat(dur, inner=[length(t)])
    offsets = repeat(vcat(0, cumsum(dur)[1:end-1]), inner=[length(t)])
    tr     .= (tr .* scale) .+ offsets
    return tr
end
function times(t, dur::Real)
    return dur .* t
end

""" unit_time """
function unit_time(tq, t, t_unit, periodic, dur::Real)
    return interpolate_times(t .* dur, t_unit, periodic, tq)
end
function unit_time(tq, t, t_unit, periodic, dur)
    return interpolate_times(times(t, dur), repeat(t_unit, length(dur)), periodic, tq)
end
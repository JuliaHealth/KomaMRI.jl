abstract type AbstractTimeSpan{T<:Real} end

"""
    timerange = TimeRange(t_start, t_end)

TimeRange struct. It is a specialized type that inherits from AbstractTimeSpan and 
defines a time interval, with start and end times.

# Arguments
- `t_start`: (`::Real`, `[s]`) start time
- `t_end`: (`::Real`, `[s]`) end time

# Returns
- `timerange`: (`::TimeRange`) TimeRange struct

# Examples
```julia-repl
julia> timerange = TimeRange(0.0, 1.0)
```
"""
@with_kw struct TimeRange{T<:Real} <: AbstractTimeSpan{T} 
   t_start  ::T 
   t_end    ::T
   @assert t_end >= t_start "t_end must be greater or equal than t_start"
end

""" Constructors """
TimeRange(t_start) = TimeRange(t_start, t_start)

""" times """
times(ts::TimeRange) = [ts.t_start, ts.t_end]

"""
    t_unit = unit_time(t, time_range)

The `unit_time` function normalizes a given array of time values t 
to a unit interval [0, 1] based on a specified start time `t_start` and end time `t_end`. 
This function is used for non-periodic motions, where each element of t is transformed 
to fit within the range [0, 1] based on the provided start and end times.

![Unit Time](../assets/unit-time.svg)

# Arguments
- `t`: (`::AbstractArray{T<:Real}`, `[s]`) array of time values to be normalized
- `time_range`: (`::TimeRange{T<:Real}`, `[s]`) time interval (defined by `t_start` and `t_end`) over which we want to normalise

# Returns
- `t_unit`: (`::AbstractArray{T<:Real}`, `[s]`) array of normalized time values

# Examples
```julia-repl
julia> t_unit = KomaMRIBase.unit_time([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], TimeRange(1.0, 4.0))
6-element Vector{Float64}:
 0.0
 0.0
 0.333
 0.666
 1.0
 1.0
```
"""
function unit_time(t, ts::TimeRange{T}) where {T<:Real}
    if ts.t_start == ts.t_end
        return (t .>= ts.t_start) .* oneunit(T)
    else
        tmp = max.((t .- ts.t_start) ./ (ts.t_end - ts.t_start), zero(T))
        return min.(tmp, oneunit(T))
    end
end


"""
    periodic = Periodic(period, asymmetry)

Periodic struct. It is a specialized type that inherits from AbstractTimeSpan, 
designed to work with time intervals that repeat periodically. It includes a measure of
asymmetry in order to recreate a asymmetric period.

# Arguments
- `period`: (`::Real`, `[s]`) period duration
- `asymmetry`: (`::Real`, `=0.5`) temporal asymmetry factor. Between 0 and 1.

# Returns
- `periodic`: (`::Periodic`) Periodic struct

# Examples
```julia-repl
julia> periodic = Periodic(1.0, 0.2)
```
"""
@with_kw struct Periodic{T<:Real} <: AbstractTimeSpan{T}
   period::T
   asymmetry::T = typeof(period)(0.5)
   t_start::T = zero(typeof(period))
end

""" Constructors """
Periodic(period) = Periodic(period, typeof(period)(0.5))

""" times """
times(ts::Periodic{T}) where {T<:Real} = ts.t_start .+ [0, ts.period * ts.asymmetry, ts.period] 

"""
    t_unit = unit_time(t, periodic)

The `unit_time` function normalizes a given array 
of time values t to a unit interval [0, 1] for periodic motions, 
based on a specified period and an asymmetry factor. 
This function is useful for creating triangular waveforms 
or normalizing time values in periodic processes.

![Unit Time Triangular](../assets/unit-time-triangular.svg)

# Arguments
- `t`: (`::AbstractArray{T<:Real}`, `[s]`) array of time values to be normalized
- `periodic`: (`::Periodic{T<:Real}`, `[s]`) information about the `period` and the temporal `asymmetry`
# Returns
- `t_unit`: (`::AbstractArray{T<:Real}`, `[s]`) array of normalized time values

# Examples
```julia-repl
julia> t_unit = KomaMRIBase.unit_time([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], Periodic(4.0, 0.5))
6-element Vector{Float64}:
 0.0
 0.5
 1.0
 0.5
 0.0
 0.5
```
"""
function unit_time(t, ts::Periodic{T}) where {T<:Real}
    t_rise = ts.period * ts.asymmetry
    t_fall = ts.period * (oneunit(T) - ts.asymmetry)
    t_relative = mod.(t, ts.period)
    if t_rise == 0
        t_unit = ifelse.(t_relative .< t_rise, zero(T), oneunit(T) .- t_relative ./ t_fall)
    elseif t_fall == 0
        t_unit = ifelse.(t_relative .< t_rise, t_relative ./ t_rise, oneunit(T))
    else
        t_unit = ifelse.( t_relative .< t_rise, t_relative ./ t_rise, oneunit(T) .- (t_relative .- t_rise) ./ t_fall)
    end
    return t_unit
end


"""
    timecurve = TimeCurve(x, y, periodic, delay)

TimeCurve struct. It is a specialized type that inherits from AbstractTimeSpan and 
defines a time curve. It is defined by the `x` and `y` vectors, an optional `periodic`
flag, and an optional delay.

# Arguments
- `x`: (`::AbstractVector{<:Real}`, `[s]`) time vector
- `y`: (`::AbstractVector{<:Real}`) y vector, it needs to be scaled between 0 and 1
- `periodic`: (`::Bool`, `=false`) indicates whether the time curve should be periodically repeated
- `delyay`: (`<:Real`, `=0.0`) initial delay which takes effect before the time curve

# Returns
- `timecurve`: (`::TimeCurve`) TimeCurve struct

# Examples
```julia-repl
julia> timecurve = TimeCurve(x=[0.0, 0.1, 0.3, 0.4], y=[0.0, 0.6, 0.2, 0.0], periodic=true, delay=0.0)
```
"""
@with_kw struct TimeCurve{T<:Real} <: AbstractTimeSpan{T} 
    x::AbstractVector{T}
    y::AbstractVector{T}
    periodic::Bool = false
    t_start::T=x[1]
end

""" times """
times(tc::TimeCurve) = tc.x
""" unit_time """
function unit_time(t, tc::TimeCurve{T}) where T<:Real
    itp = GriddedInterpolation((tc.x, ), tc.y, Gridded(Linear()))
    return extrapolate(itp, tc.periodic ? Interpolations.Periodic() : Flat()).(t)
end
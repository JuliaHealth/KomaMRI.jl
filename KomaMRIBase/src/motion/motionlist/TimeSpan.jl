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
end

""" Constructors """
Periodic(period) = Periodic(period, typeof(period)(0.5))

""" times """
times(ts::Periodic{T}) where {T<:Real} = [zero(T), ts.period * ts.asymmetry, ts.period] 
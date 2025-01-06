"""
    timecurve = TimeCurve(t, t_unit, periodic, periods)

TimeCurve struct. It is a specialized type that defines a time curve, which represents 
the temporal behavior of motion. This curve is defined by two vectors: 
`t` and `t_unit`, which represent the horizontal (x-axis) and vertical (y-axis) axes 
of the curve, respectively. To some extent, this curve can be associated with animation curves,
commonly found in software for video editing, 3D scene creation, or video game development.

Additionally, the TimeCurve struct contains two more fields, independent of each other:
`periodic` is a Boolean that indicates whether the time curve should be repeated periodically.
`periods` contains as many elements as repetitions are desired in the time curve. 
Each element specifies the scaling factor for that repetition.

# Arguments
- `t`: (`::AbstractVector{<:Real}`, `[s]`) time vector
- `t_unit`: (`::AbstractVector{<:Real}`) y vector, it needs to be scaled between 0 and 1
- `periodic`: (`::Bool`, `=false`) indicates whether the time curve should be periodically repeated
- `periods`: (`::Union{<:Real,AbstractVector{<:Real}}`, `=1.0`): represents the relative duration 
    of each period with respect to the baseline duration defined by `t[end] - t[1]`. 
    In other words, it acts as a scaling factor to lengthen or shorten specific periods. 
    This allows for the creation of patterns such as arrhythmias or other variations in periodicity.

# Returns
- `timecurve`: (`::TimeCurve`) TimeCurve struct

# Examples
1. Non-periodic motion with a single repetition:
    ```julia-repl
    julia> timecurve = TimeCurve(t=[0.0, 0.2, 0.4, 0.6], t_unit=[0.0, 0.2, 0.5, 1.0])
    ```
    ![Time Curve 1](../assets/time-curve-1.svg)
2. Periodic motion with a single repetition:
    ```julia-repl
    julia> timecurve = TimeCurve(t=[0.0, 0.2, 0.4, 0.6], t_unit=[0.0, 1.0, 1.0, 0.0], periodic=true)
    ```
    ![Time Curve 2](../assets/time-curve-2.svg)
3. Non-periodic motion with multiple repetitions:
    ```julia-repl
    julia> timecurve = TimeCurve(t=[0.0, 0.2, 0.4, 0.6], t_unit=[0.0, 1.0, 1.0, 0.0], periods=[1.0, 0.5, 1.5])
    ```
    ![Time Curve 3](../assets/time-curve-3.svg)
4. Periodic motion with multiple repetitions:
    ```julia-repl
    julia> timecurve = TimeCurve(t=[0.0, 0.2, 0.4, 0.6], t_unit=[0.0, 1.0, 1.0, 0.0], periods=[1.0, 0.5, 1.5], periodic=true)
    ```
    ![Time Curve 4](../assets/time-curve-4.svg)
"""
@with_kw struct TimeCurve{T<:Real}
    t::AbstractVector{T}
    t_unit::AbstractVector{T}
    periodic::Bool                       = false
    periods::Union{T,AbstractVector{T}} = oneunit(eltype(t))
    t_start::T                           = t[1]
    t_end::T                             = t[end]
    @assert check_unique(t) "Vector t=$(t) contains duplicate elements. Please ensure all elements in t are unique and try again"
end

check_unique(t) = true
check_unique(t::Vector) = length(t) == length(unique(t))

# Main Constructors
TimeCurve(t, t_unit, periodic, periods) = TimeCurve(t=t, t_unit=t_unit, periodic=periodic, periods=periods)
TimeCurve(t, t_unit) = TimeCurve(t=t, t_unit=t_unit)
# Custom constructors
"""
    timerange = TimeRange(t_start, t_end)

The `TimeRange` function is a custom constructor for the `TimeCurve` struct. 
It allows defining a simple time interval, with start and end times.

# Arguments
- `t_start`: (`::Real`, `[s]`, `=0.0`) start time
- `t_end`: (`::Real`, `[s]`, `=1.0`) end time

# Returns
- `timerange`: (`::TimeCurve`) TimeCurve struct

# Examples
```julia-repl
julia> timerange = TimeRange(t_start=0.6, t_end=1.4)
```
![Time Range](../assets/time-range.svg)
"""
TimeRange(t_start::T, t_end::T) where T = TimeCurve(t=[t_start, t_end], t_unit=[zero(T), oneunit(T)])
TimeRange(; t_start=0.0, t_end=1.0)     = TimeRange(t_start, t_end)
"""
    periodic = Periodic(period, asymmetry)

The `Periodic` function is a custom constructor for the `TimeCurve` struct.
It allows defining time intervals that repeat periodically with a triangular period. 
It includes a measure of asymmetry in order to recreate a asymmetric period.

# Arguments
- `period`: (`::Real`, `[s]`, `=1.0`) period duration
- `asymmetry`: (`::Real`, `=0.5`) temporal asymmetry factor. Between 0 and 1.

# Returns
- `periodic`: (`::TimeCurve`) TimeCurve struct

# Examples
```julia-repl
julia> periodic = Periodic(period=1.0, asymmetry=0.2)
```
![Periodic](../assets/periodic.svg)
"""
Periodic(period::T, asymmetry::T) where T = TimeCurve(t=[zero(T), period*asymmetry, period], t_unit=[zero(T), oneunit(T), zero(T)])
Periodic(; period=1.0, asymmetry=0.5)     = Periodic(period, asymmetry)

""" Compare two TimeCurves """
Base.:(==)(t1::TimeCurve, t2::TimeCurve) = reduce(&, [getfield(t1, field) == getfield(t2, field) for field in fieldnames(typeof(t1))])
Base.:(≈)(t1::TimeCurve, t2::TimeCurve)  = reduce(&, [getfield(t1, field)  ≈ getfield(t2, field) for field in fieldnames(typeof(t1))])

""" times & unit_time """
# Although the implementation of these two functions when `per` is a vector is valid 
# for all cases, it performs unnecessary and costly operations when `per` is a scalar. 
# Therefore, it has been decided to use method dispatch between these two cases.
function times(t, per::Real)
    return per .* t
end
function times(t, per::AbstractVector)
    tr      = repeat(t, length(per))
    scale   = repeat(per, inner=[length(t)])
    offsets = repeat(vcat(0, cumsum(per)[1:end-1]), inner=[length(t)])
    tr     .= (tr .* scale) .+ offsets
    return tr
end
function unit_time(tq, t, t_unit, periodic, per::Real)
    return interpolate_times(t .* per, t_unit, periodic, tq)
end
function unit_time(tq, t, t_unit, periodic, per::AbstractVector)
    return interpolate_times(times(t, per), repeat(t_unit, length(per)), periodic, tq)
end
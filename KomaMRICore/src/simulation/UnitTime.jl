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
function unit_time(t::AbstractArray{T}, ts::TimeRange{T}) where {T<:Real}
    if ts.t_start == ts.t_end
        return (t .>= ts.t_start) .* oneunit(T)
    else
        tmp = max.((t .- ts.t_start) ./ (ts.t_end - ts.t_start), zero(T))
        t = min.(tmp, oneunit(T))
        # _ = sum(t) # Dummy (oneAPI bug)
        KA.synchronize(KA.get_backend(t))
        return t
    end
end
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
function unit_time(t::AbstractArray{T}, ts::Periodic{T}) where {T<:Real}
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

"""
    t_unit = unit_time(t, t_start, t_end)

The `unit_time` function normalizes a given array of time values t 
to a unit interval [0, 1] based on a specified start time t_start and end time t_end. 
This function is used for non-periodic motions, where each element of t is transformed 
to fit within the range [0, 1] based on the provided start and end times.

![Unit Time](../assets/unit-time.svg)

# Arguments
- `t`: (`::AbstractArray{T<:Real}`, `[s]`) array of time values to be normalized
- `t_start`: (`::T`, `[s]`) start time for normalization
- `t_end`: (`::T`, `[s]`) end time for normalization

# Returns
- `t_unit`: (`::AbstractArray{T<:Real}`, `[s]`) array of normalized time values

# Examples
```julia-repl
julia> t_unit = KomaMRIBase.unit_time([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 1.0, 4.0)
6-element Vector{Float64}:
 0.0
 0.0
 0.3333333333333333
 0.6666666666666666
 1.0
 1.0
"""
function unit_time(t::AbstractArray{T}, t_start::T, t_end::T) where {T<:Real}
    if t_start == t_end
        return (t .>= t_start) .* oneunit(T)
    else
        return min.(max.((t .- t_start) ./ (t_end - t_start), zero(T)), oneunit(T))
    end
end

"""
    t_unit = unit_time_triangular(t, period, asymmetry)

The `unit_time_triangular` function normalizes a given array 
of time values t to a unit interval [0, 1] for periodic motions, 
based on a specified period and an asymmetry factor. 
This function is useful for creating triangular waveforms 
or normalizing time values in periodic processes.

![Unit Time Triangular](../assets/unit-time-triangular.svg)

# Arguments
- `t`: (`::AbstractArray{T<:Real}`, `[s]`) array of time values to be normalized
- `period`: (`::T`, `[s]`) the period of the triangular waveform
- `asymmetry`: (`::T`) asymmetry factor, a value in the range (0, 1) indicating the fraction of the period in the rising part of the triangular wave

# Returns
- `t_unit`: (`::AbstractArray{T<:Real}`, `[s]`) array of normalized time values

# Examples
```julia-repl
julia> t_unit = KomaMRIBase.unit_time_triangular([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 4.0, 0.5)
6-element Vector{Float64}:
 0.0
 0.5
 1.0
 0.5
 0.0
 0.5
"""
function unit_time_triangular(t::AbstractArray{T}, period::T, asymmetry::T) where {T<:Real}
    t_rise = period * asymmetry
    t_fall = period * (1 - asymmetry)
    t_relative = mod.(t, period)
    t_unit =
        ifelse.(
            t_relative .< t_rise,
            t_relative ./ t_rise,
            1 .- (t_relative .- t_rise) ./ t_fall,
        )
    return t_unit
end

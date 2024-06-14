"""
    t_unit = unit_time(t, t_start, t_end)

For non-periodic motions
"""
function unit_time(t::AbstractArray{T}, t_start::T, t_end::T) where {T<:Real}
    if t_start == t_end
        return (t .>= t_start) .* one(T)
    else
        return min.(max.((t .- t_start) ./ (t_end - t_start), zero(T)), one(T))
    end
end

"""
    t_unit = unit_time_triangular(t, period, asymmetry)

For periodic motions
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

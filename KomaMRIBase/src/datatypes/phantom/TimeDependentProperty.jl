struct TimeDependentProperty{T<:Real}
    value::AbstractMatrix{T}
    time::TimeCurve{T}
end

get_property_value(p::TimeDependentProperty) = p.value[:, 1]

"""Spin-resolved phantom properties at simulation times `t` (vector or matrix columns)."""
function get_spin_properties(p, t)
    return (
        ρ=get_spin_property(p.ρ, t),
        T1=get_spin_property(p.T1, t),
        T2=get_spin_property(p.T2, t),
        T2s=get_spin_property(p.T2s, t),
        Δw=get_spin_property(p.Δw, t),
        Dλ1=get_spin_property(p.Dλ1, t),
        Dλ2=get_spin_property(p.Dλ2, t),
        Dθ=get_spin_property(p.Dθ, t),
    )
end

"""Block-end spin properties used for GPU kernels and block-final relaxation."""
function get_spin_properties_block_end(p, t)
    props = get_spin_properties(p, t)
    return (
        ρ=get_spin_property_at_end(props.ρ),
        T1=get_spin_property_at_end(props.T1),
        T2=get_spin_property_at_end(props.T2),
        T2s=get_spin_property_at_end(props.T2s),
        Δw=get_spin_property_at_end(props.Δw),
        Dλ1=get_spin_property_at_end(props.Dλ1),
        Dλ2=get_spin_property_at_end(props.Dλ2),
        Dθ=get_spin_property_at_end(props.Dθ),
    )
end

get_spin_property(p::AbstractVector, t) = p
function get_spin_property(p::TimeDependentProperty, t)
    Ns = size(p.value, 1)
    t_unit = unit_time(t, p.time)
    return resample(interpolate(p.value, Gridded(Linear()), Val(Ns), t_unit), t_unit)
end

"""
    get_spin_property_at_end(p::AbstractMatrix)

Returns spin-wise property values at the last queried time in `p`, with shape
`(Nspins, Ntimes)`.
"""
get_spin_property_at_end(p::AbstractMatrix) = vec(selectdim(p, 2, size(p, 2)))
get_spin_property_at_end(p::AbstractVector) = p

Base.:(==)(p1::TimeDependentProperty, p2::TimeDependentProperty) = p1.value == p2.value && p1.time == p2.time
Base.:(≈)(p1::TimeDependentProperty,  p2::TimeDependentProperty) = p1.value  ≈ p2.value && p1.time  ≈ p2.time
Base.:(==)(::TimeDependentProperty, ::AbstractVector) = false
Base.:(==)(::AbstractVector, ::TimeDependentProperty) = false
Base.:(≈)(::TimeDependentProperty, ::AbstractVector)  = false
Base.:(≈)(::AbstractVector, ::TimeDependentProperty)  = false

Base.getindex(p::TimeDependentProperty, i) = TimeDependentProperty(p.value[i, :], p.time)
Base.view(p::TimeDependentProperty, i) = TimeDependentProperty(@view(p.value[i, :]), p.time)

times(::AbstractVector) = Float64[]
times(p::TimeDependentProperty) = times(p.time)

*(α::Real, p::TimeDependentProperty) = TimeDependentProperty(α .* p.value, p.time)

vcat_property(p1::AbstractVector, p2::AbstractVector) = [p1; p2]
vcat_property(p1::TimeDependentProperty, p2::AbstractVector) = TimeDependentProperty([p1.value; repeat(p2, 1, size(p1.value, 2))], p1.time)
vcat_property(p1::AbstractVector, p2::TimeDependentProperty) = TimeDependentProperty([repeat(p1, 1, size(p2.value, 2)); p2.value], p2.time)
vcat_property(p1::TimeDependentProperty, p2::TimeDependentProperty) = begin
    @assert p1.time == p2.time "Time ranges must be the same"
    return TimeDependentProperty([p1.value; p2.value], p1.time)
end
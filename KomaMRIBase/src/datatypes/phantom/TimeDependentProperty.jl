struct TimeDependentProperty{T<:Real}
    value::AbstractMatrix{T}
    time::TimeCurve{T}
end

get_property_value(p::TimeDependentProperty) = p.value[:, 1]

get_spin_property(p::AbstractVector, t) = p
function get_spin_property(p::TimeDependentProperty, t)
    Ns = size(p.value, 1)
    t_unit = unit_time(t, p.time)
    return resample(interpolate(p.value, Gridded(Linear()), Val(Ns), t_unit), t_unit)
end

get_spin_property_at_end(p::AbstractMatrix) = selectdim(p, 2, size(p, 2))
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
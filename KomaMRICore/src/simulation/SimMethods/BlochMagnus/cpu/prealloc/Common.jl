abstract type BlochMagnusCPUPrealloc <: PreallocResult end

cbuf(obj::Phantom) = zero.(complex.(obj.ρ))
rbuf(obj::Phantom) = zero.(obj.ρ)
off_resonance_buffer(obj::Phantom) = obj.Δw ./ eltype(obj.ρ)(2π .* γ)

function Base.view(p::P, i::UnitRange) where {P<:BlochMagnusCPUPrealloc}
    fields = ntuple(j -> prealloc_view(getfield(p, j), i), Val(fieldcount(P)))
    return Base.typename(P).wrapper(fields...)
end

prealloc_view(x::AbstractVector, i) = view(x, i)
prealloc_view(x::AbstractMatrix, i) = view(x, i, :)
prealloc_view(coordinates::MotionCoordinates, i) =
    view_motion_coordinates(coordinates, i)
prealloc_view(sensitivities::CoilSensitivities, i) = CoilSensitivities(
    @view(sensitivities.values[i, :]),
    sensitivities.interpolators,
)
prealloc_view(x, _) = x

precession_buffers(p::BlochMagnusCPUPrealloc) = p.ωz_0, p.ωz_1

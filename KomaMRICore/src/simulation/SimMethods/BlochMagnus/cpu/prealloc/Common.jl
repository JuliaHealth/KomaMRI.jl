abstract type BlochMagnusCPUPrealloc{T} <: PreallocResult{T} end

cbuf(obj::Phantom{T}) where {T<:Real} = zeros(Complex{T}, size(obj.x))
rbuf(obj::Phantom{T}) where {T<:Real} = zeros(T, size(obj.x))
off_resonance_buffer(obj::Phantom{T}) where {T<:Real} = obj.Δw ./ T(2π .* γ)

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

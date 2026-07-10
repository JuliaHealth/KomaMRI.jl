abstract type BlochMagnusCPUPrealloc <: PreallocResult end

cbuf(obj::Phantom) = zeros(Complex{eltype(obj.ρ)}, size(obj.x))
rbuf(obj::Phantom) = zeros(eltype(obj.ρ), size(obj.x))
off_resonance_buffer(obj::Phantom) = obj.Δw ./ eltype(obj.ρ)(2π .* γ)

function Base.view(p::P, i::UnitRange) where {P<:BlochMagnusCPUPrealloc}
    fields = ntuple(j -> view(getfield(p, j), i), Val(fieldcount(P)))
    return Base.typename(P).wrapper(fields...)
end

precession_buffers(p::BlochMagnusCPUPrealloc) = p.ωz_0, p.ωz_1

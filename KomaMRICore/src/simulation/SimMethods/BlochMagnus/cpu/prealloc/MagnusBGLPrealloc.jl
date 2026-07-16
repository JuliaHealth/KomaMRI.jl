struct BlochMagnusBGLCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T}
} <: BlochMagnusCPUPrealloc{T}
    ωxy_minus::CV
    ωz_minus::RV
    ωxy_center::CV
    ωz_center::RV
    ωxy_plus::CV
    ωz_plus::RV
    i0xy::CV
    i0z::RV
    i1xy::CV
    i1z::RV
    i2xy::CV
    i2z::RV
    jxy::CV
    jz::RV
    boxxy::CV
    boxz::RV
    θxy::CV
    θz::RV
    rotation_norm::RV
    α::CV
    β::CV
    ΔBz::RV
    Maux_xy::CV
    Maux_z::RV
end

prealloc(sim_method::BlochMagnusBGL4, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    BlochMagnusBGLCPUPrealloc(
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
    )

prealloc(sim_method::BlochMagnusBGL6, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    prealloc(BlochMagnusBGL4(), backend, obj, M, max_block_length, groupsize)

precession_buffers(p::BlochMagnusBGLCPUPrealloc) = p.ωz_minus, p.ωz_plus

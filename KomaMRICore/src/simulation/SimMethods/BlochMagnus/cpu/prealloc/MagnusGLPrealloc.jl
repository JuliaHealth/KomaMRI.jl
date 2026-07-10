struct BlochMagnusGLCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T}
} <: BlochMagnusCPUPrealloc{T}
    ωxy_minus::CV
    ωz_minus::RV
    ωxy_plus::CV
    ωz_plus::RV
    θxy::CV
    θz::RV
    rotation_norm::RV
    α::CV
    β::CV
    ΔBz::RV
    Maux_xy::CV
    Maux_z::RV
    relaxation::RelaxationCPUPrealloc{T,RV}
end

prealloc(sim_method::BlochMagnusGL4, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    BlochMagnusGLCPUPrealloc(
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
        relaxation_prealloc(obj),
    )

prealloc(sim_method::BlochMagnusGL2, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    prealloc(BlochMagnusGL4(), backend, obj, M, max_block_length, groupsize)

precession_buffers(p::BlochMagnusGLCPUPrealloc) = p.ωz_minus, p.ωz_plus

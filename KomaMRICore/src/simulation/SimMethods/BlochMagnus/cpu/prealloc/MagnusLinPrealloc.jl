struct BlochMagnusLinCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T}
} <: BlochMagnusCPUPrealloc{T}
    ωxy_0::CV
    ωz_0::RV
    ωxy_1::CV
    ωz_1::RV
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

prealloc(sim_method::BlochMagnusLin2, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    BlochMagnusLinCPUPrealloc(
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
        relaxation_prealloc(obj),
    )

prealloc(sim_method::BlochMagnusLinComm2, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    prealloc(BlochMagnusLin2(), backend, obj, M, max_block_length, groupsize)

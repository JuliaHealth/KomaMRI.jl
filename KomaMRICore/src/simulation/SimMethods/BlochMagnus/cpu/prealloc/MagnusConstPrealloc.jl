struct BlochMagnusConstCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T}
} <: BlochMagnusCPUPrealloc{T}
    ωxy_0::CV
    ωz_0::RV
    ωz_1::RV
    θxy::CV
    θz::RV
    rotation_norm::RV
    α::CV
    β::CV
    ΔBz::RV
    Maux_xy::CV
    Maux_z::RV
end

prealloc(sim_method::BlochMagnusConst1, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize, sys::Scanner) where {T<:Real} =
    BlochMagnusConstCPUPrealloc(
        cbuf(obj), rbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
    )

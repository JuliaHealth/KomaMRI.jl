struct BlochMagnusMidCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T}
} <: BlochMagnusCPUPrealloc{T}
    ωxy_m::CV
    ωz_m::RV
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

prealloc(sim_method::BlochMagnusMid2, backend::KA.CPU, obj::Phantom{T}, M::Mag{T}, max_block_length::Integer, groupsize) where {T<:Real} =
    BlochMagnusMidCPUPrealloc(
        cbuf(obj), rbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
    )

precession_buffers(p::BlochMagnusMidCPUPrealloc) = p.ωz_m, p.ωz_1

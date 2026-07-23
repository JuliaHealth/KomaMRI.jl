struct BlochMagnusQuadCPUPrealloc{CV<:AbstractVector,RV<:AbstractVector} <: BlochMagnusCPUPrealloc
    ωxy_0::CV
    ωz_0::RV
    ωxy_m::CV
    ωz_m::RV
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
end

prealloc(sim_method::BlochMagnusQuad4, backend::KA.CPU, obj::Phantom, M::Mag, max_block_length::Integer, groupsize) =
    BlochMagnusQuadCPUPrealloc(
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
    )

prealloc(sim_method::BlochMagnusQuad2, backend::KA.CPU, obj::Phantom, M::Mag, max_block_length::Integer, groupsize) =
    prealloc(BlochMagnusQuad4(), backend, obj, M, max_block_length, groupsize)

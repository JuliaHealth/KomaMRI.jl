struct BlochMagnusMidCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T},S,P
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
    sens::S
    coordinates::P
end

prealloc(::BlochMagnusMid2, backend::KA.CPU, obj, M, max_block_length, _max_adc_samples, _groupsize, sys) =
    BlochMagnusMidCPUPrealloc(
        cbuf(obj), rbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
        prealloc_sensitivities(sys.receiver, obj),
        prealloc_motion_coordinates(obj.motion, backend, obj, max_block_length),
    )

precession_buffers(p::BlochMagnusMidCPUPrealloc) = p.ωz_m, p.ωz_1

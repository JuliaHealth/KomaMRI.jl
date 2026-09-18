struct BlochMagnusConstCPUPrealloc{
    T,CV<:AbstractVector{Complex{T}},RV<:AbstractVector{T},S,P,R
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
    sens::S
    coordinates::P
    spin_reset::R
end

prealloc(::BlochMagnusConst1, backend::KA.CPU, obj, M, max_block_length, _max_adc_samples, _groupsize, sys) =
    BlochMagnusConstCPUPrealloc(
        cbuf(obj), rbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
        prealloc_sensitivities(sys.receiver, obj),
        prealloc_motion_coordinates(obj.motion, backend, obj, max_block_length),
        prealloc_spin_reset(obj.motion, backend, length(obj), max_block_length),
    )

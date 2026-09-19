struct BlochMagnusGLCPUPrealloc{
    CV<:AbstractVector,RV<:AbstractVector,S,P
} <: BlochMagnusCPUPrealloc
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
    sens::S
    coordinates::P
end

prealloc(::BlochMagnusGL4, backend::KA.CPU, obj, M, max_block_length, _max_adc_samples, _groupsize, sys) =
    BlochMagnusGLCPUPrealloc(
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
        prealloc_sensitivities(sys.receiver, obj),
        prealloc_motion_coordinates(obj.motion, backend, obj, max_block_length),
    )

prealloc(::BlochMagnusGL2, backend::KA.CPU, obj, M, max_block_length, max_adc_samples, groupsize, sys) =
    prealloc(BlochMagnusGL4(), backend, obj, M, max_block_length, max_adc_samples, groupsize, sys)

precession_buffers(p::BlochMagnusGLCPUPrealloc) = p.ωz_minus, p.ωz_plus

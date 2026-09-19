struct BlochMagnusLinCPUPrealloc{
    CV<:AbstractVector,RV<:AbstractVector,S,P
} <: BlochMagnusCPUPrealloc
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
    sens::S
    coordinates::P
end

prealloc(::BlochMagnusLin2, backend::KA.CPU, obj, M, max_block_length, _max_adc_samples, _groupsize, sys) =
    BlochMagnusLinCPUPrealloc(
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj),
        cbuf(obj), rbuf(obj), rbuf(obj),
        similar(M.xy), similar(M.xy),
        off_resonance_buffer(obj),
        similar(M.xy), similar(M.z),
        prealloc_sensitivities(sys.receiver, obj),
        prealloc_motion_coordinates(obj.motion, backend, obj, max_block_length),
    )

prealloc(::BlochMagnusLinComm2, backend::KA.CPU, obj, M, max_block_length, max_adc_samples, groupsize, sys) =
    prealloc(BlochMagnusLin2(), backend, obj, M, max_block_length, max_adc_samples, groupsize, sys)

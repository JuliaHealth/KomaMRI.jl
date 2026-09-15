function combine_coil_signal!(sig, Mxy::AbstractMatrix, sens, ::KA.CPU)
    sig .= transpose(Mxy) * sens
    return nothing
end

function combine_coil_signal!(sig, Mxy::AbstractVector, sens, ::KA.CPU)
    sig .= vec(transpose(Mxy) * sens)
    return nothing
end

function combine_coil_signal!(sig, Mxy::AbstractVector, sens, ::KA.GPU)
    sig .= vec(sum(reshape(Mxy, :, 1) .* sens; dims=1))
    return nothing
end

function combine_coil_signal!(sig, Mxy::AbstractMatrix, sens, backend::KA.GPU)
    for sample in axes(Mxy, 2)
        combine_coil_signal!(
            @view(sig[sample, :]), @view(Mxy[:, sample]), sens, backend,
        )
    end
    return nothing
end

combine_coil_signal!(sig, Mxy, sens) =
    combine_coil_signal!(sig, Mxy, sens, KA.get_backend(sig))

KomaMRIBase.get_sens(receiver, x, y, z, ::KA.Backend) =
    get_sens(receiver, x, y, z)

struct CoilSensitivities{V,I}
    values::V
    interpolators::I
end

prealloc_sensitivities(receiver::UniformCoilSens, _) = receiver
prealloc_sensitivities(receiver::AbstractRFReceiveSystem, obj) =
    CoilSensitivities(get_sens(receiver, obj.x, obj.y, obj.z), nothing)
prealloc_sensitivities(receiver::ArbitraryCoilSens, obj) = CoilSensitivities(
    get_sens(receiver, obj.x, obj.y, obj.z),
    KomaMRIBase.coil_interpolators(receiver),
)

update_sensitivities!(::UniformCoilSens, _, _, _) = nothing
update_sensitivities!(::CoilSensitivities, _, _, ::NoMotion) = nothing
function update_sensitivities!(
    sensitivities::CoilSensitivities, receiver::AbstractRFReceiveSystem,
    (x, y, z), ::Union{Motion,MotionList},
)
    KomaMRIBase.get_sens!(sensitivities.values, receiver, x, y, z)
    return nothing
end
function update_sensitivities!(
    sensitivities::CoilSensitivities, receiver::ArbitraryCoilSens,
    (x, y, z), ::Union{Motion,MotionList},
)
    KomaMRIBase.get_sens!(
        sensitivities.values, receiver, x, y, z, sensitivities.interpolators,
    )
    return nothing
end

function acquire_signal!(sig, Mxy::AbstractVector, ::UniformCoilSens, _)
    sig .= sum(Mxy)
    return nothing
end

function acquire_signal!(sig, Mxy::AbstractMatrix, ::UniformCoilSens, _)
    sig .= vec(sum(Mxy; dims=1))
    return nothing
end

function acquire_signal!(sig, Mxy, sensitivities::CoilSensitivities, _)
    combine_coil_signal!(sig, Mxy, sensitivities.values)
    return nothing
end

function acquire_signal!(sig, Mxy, receiver, (x, y, z))
    combine_coil_signal!(
        sig, Mxy, get_sens(receiver, x, y, z, KA.get_backend(sig)),
    )
    return nothing
end

# Uniform Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, _, ::UniformCoilSens, Mxy::AbstractMatrix, ::NoMotion, _, _)
    sig .= @views transpose(sum(Mxy; dims=1))
    return nothing
end

function acquire_signal!(
    sig, _, ::UniformCoilSens, Mxy::AbstractMatrix,
    ::Union{Motion,MotionList}, _, _,
)
    sig .= @views transpose(sum(Mxy; dims=1))
    return nothing
end

# Uniform other sim methods
function acquire_signal!(
    sig, _phantom, ::UniformCoilSens, Mxy::AbstractVector,
    ::NoMotion, _coords, _sens=nothing,
)
    sig[1] = sum(Mxy)
    return nothing
end

function acquire_signal!(
    sig, _phantom, ::UniformCoilSens, Mxy::AbstractVector,
    ::Union{Motion,MotionList}, _coords, _sens=nothing,
)
    sig[1] = sum(Mxy)
    return nothing
end

# Nonuniform Bloch Simple run_spin_precession not sample based, without motion
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
    combine_coil_signal!(sig, Mxy, sens, KA.get_backend(Mxy))

function acquire_signal!(sig, p, receiver::AbstractRFReceiveSystem, Mxy::AbstractMatrix, ::NoMotion, _, _)
    sens = receiver_sensitivities(receiver, p.x, p.y, p.z)
    combine_coil_signal!(sig, Mxy, sens)
    return nothing
end

# Nonuniform Bloch Simple run_spin_precession not sample based, with motion
function acquire_signal!(
    sig, p, receiver::AbstractRFReceiveSystem, Mxy::AbstractMatrix,
    ::Union{Motion,MotionList}, (x, y, z), adc,
)
    for (sample, time_index) in enumerate(adc .+ 1)
        coords = (@view(x[:, time_index]), @view(y[:, time_index]), @view(z[:, time_index]))
        acquire_signal!(
            @view(sig[sample, :]), p, receiver, Mxy[:, sample],
            p.motion, coords,
        )
    end
    return nothing
end

# Nonuniform other sim methods without motion
function acquire_signal!(
    sig, p, receiver::AbstractRFReceiveSystem, Mxy::AbstractVector,
    ::NoMotion, _coords, sens=receiver_sensitivities(receiver, p.x, p.y, p.z),
)
    combine_coil_signal!(sig, Mxy, sens)
    return nothing
end

function acquire_signal!(
    sig, _phantom, receiver::AbstractRFReceiveSystem, Mxy::AbstractVector,
    ::Union{Motion,MotionList}, coords, _sens=nothing,
)
    combine_coil_signal!(sig, Mxy, receiver_sensitivities(receiver, coords...))
    return nothing
end

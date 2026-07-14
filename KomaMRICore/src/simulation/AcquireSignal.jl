# Uniform Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, _, ::UniformCoilSens, Mxy::AbstractMatrix)
    sig .= @views transpose(sum(Mxy; dims=1))
    return nothing
end

# Uniform other sim methods
function acquire_signal!(sig, _, ::UniformCoilSens, Mxy::AbstractVector)
    sig .= sum(Mxy)
    return nothing
end

# Birdcage Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, p, receiver::BirdcageCoilSens, Mxy::AbstractMatrix)
    sig .= transpose(Mxy) * get_sens(receiver, p.x, p.y, p.z)
    return nothing
end

# Birdcage other sim methods
function acquire_signal!(sig, p, receiver::BirdcageCoilSens, Mxy::AbstractVector)
    sig .= vec(transpose(Mxy) * get_sens(receiver, p.x, p.y, p.z))
    return nothing
end

# Arbitrary Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, p, receiver::ArbitraryCoilSens, Mxy::AbstractMatrix)
    interpolated_coil_sens = similar(receiver.coil_sens, length(p.x), size(receiver.coil_sens, 4))
    for i in 1:size(receiver.coil_sens, 4)
        base_itp = KomaMRIBase.GriddedInterpolation((receiver.x, receiver.y, receiver.z), receiver.coil_sens[:,:,:,i], KomaMRIBase.Gridded(KomaMRIBase.Linear()))
        itp = KomaMRIBase.extrapolate(base_itp, 0f0)
        interpolated_coil_sens[:,i] = itp.(p.x, p.y, p.z)
        sig[:, i] .= vec(sum(interpolated_coil_sens[:, i] .* Mxy, dims=1))
    end
    return nothing
end

# Arbitrary other sim methods
function acquire_signal!(sig, p, receiver::ArbitraryCoilSens, Mxy::AbstractVector)
    interpolated_coil_sens = similar(receiver.coil_sens, length(p.x), size(receiver.coil_sens, 4))
    for i in 1:size(receiver.coil_sens, 4)
        base_itp = KomaMRIBase.GriddedInterpolation((receiver.x, receiver.y, receiver.z), receiver.coil_sens[:,:,:,i], KomaMRIBase.Gridded(KomaMRIBase.Linear()))
        itp = KomaMRIBase.extrapolate(base_itp, 0f0)
        interpolated_coil_sens[:,i] = itp.(p.x, p.y, p.z)
        sig[i] = sum(interpolated_coil_sens[:, i] .* Mxy)
    end
    return nothing
end

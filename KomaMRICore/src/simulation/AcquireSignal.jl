# Uniform Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, _, ::UniformCoilSens, Mxy::AbstractMatrix)
    sig .= @views transpose(sum(Mxy; dims=1))
    return nothing
end

# Uniform other sim methods
function acquire_signal!(sig, _, ::UniformCoilSens, Mxy::AbstractVector)
    sig[1] = sum(Mxy)
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
    sig .= transpose(Mxy) * get_sens(receiver, p.x, p.y, p.z)
    return nothing
end

# Arbitrary other sim methods
function acquire_signal!(sig, p, receiver::ArbitraryCoilSens, Mxy::AbstractVector)
    sig .= vec(transpose(Mxy) * get_sens(receiver, p.x, p.y, p.z))
    return nothing
end

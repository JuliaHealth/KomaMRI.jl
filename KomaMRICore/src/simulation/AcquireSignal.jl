# Uniform Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, _, ::UniformCoilSens, Mxy::AbstractMatrix, _, _, adc)
    sig .= @views transpose(sum(Mxy; dims=1))
    return nothing
end

# Uniform other sim methods
function acquire_signal!(sig, p, ::UniformCoilSens, Mxy::AbstractVector, motion, coords, sens=nothing)
    sig[1] = sum(Mxy)
    return nothing
end

# Nonuniform Bloch Simple run_spin_precession not sample based, without motion
function acquire_signal!(sig, p, receiver::Union{BirdcageCoilSens,ArbitraryCoilSens}, Mxy::AbstractMatrix, ::NoMotion, _, adc)
    sens = reshape(get_sens(receiver, p.x, p.y, p.z), size(Mxy, 1), 1, get_n_coils(receiver))
    sig .= dropdims(sum(Mxy .* sens; dims=1); dims=1)
    return nothing
end

# Nonuniform Bloch Simple run_spin_precession not sample based, with motion
function acquire_signal!(sig, p, receiver::Union{BirdcageCoilSens,ArbitraryCoilSens}, Mxy::AbstractMatrix, ::Union{Motion,MotionList}, (x, y, z), adc)
    @views sens = reshape(get_sens(receiver, vec(x[:, adc .+ 1]), vec(y[:, adc .+ 1]), vec(z[:, adc .+ 1])), size(Mxy, 1), length(adc), get_n_coils(receiver))
    sig .= dropdims(sum(Mxy .* sens; dims=1); dims=1)
    return nothing
end

# Nonuniform other sim methods without motion
function acquire_signal!(sig, p, receiver::Union{BirdcageCoilSens,ArbitraryCoilSens}, Mxy::AbstractVector, motion::NoMotion, coords, sens=get_sens(receiver, p.x, p.y, p.z))
    sig .= vec(transpose(Mxy) * sens)
    return nothing
end

# Nonuniform other sim methods with motion
function acquire_signal!(sig, p, receiver::Union{BirdcageCoilSens,ArbitraryCoilSens}, Mxy::AbstractVector, motion::Union{Motion,MotionList}, coords, sens=nothing)
    sig .= vec(transpose(Mxy) * get_sens(receiver, coords...))
    return nothing
end

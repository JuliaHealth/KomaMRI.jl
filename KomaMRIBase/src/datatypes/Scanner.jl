include("scanner/GradientCoils.jl")
include("scanner/RFCoilsTx.jl")
include("scanner/RFCoilsRx.jl")
include("scanner/HardwareLimits.jl")

export HardwareLimits, Scanner,
    AbstractGradientSystem, GradientCoils, LinearXYZ,
    AbstractReceiveSystem, RFCoilsRx, UniformCoilSens, BirdcageCoilSens,
    AbstractTransmitSystem, RFCoilsTx, UniformTransmit,
    get_n_coils, get_sens

"""
    sys = Scanner(; limits=HardwareLimits(), gradient=LinearXYZ(),
        receiver=UniformCoilSens(), transmitter=UniformTransmit())

The Scanner struct. It contains hardware limitations of the MRI resonator. It is an input
for the simulation.

# Arguments
- `limits`: (`=HardwareLimits()`) scanner hardware limits
- `gradient`: (`=LinearXYZ()`) gradient coil model
- `receiver`: (`=UniformCoilSens()`) receive coil sensitivity model
- `transmitter`: (`=UniformTransmit()`) transmit coil model

# Returns
- `sys`: (`::Scanner`) Scanner struct

# Examples
```julia-repl
julia> sys = Scanner()

julia> sys.limits.B0
```
"""
mutable struct Scanner{
    G<:AbstractGradientSystem,
    RX<:AbstractReceiveSystem,
    TX<:AbstractTransmitSystem,
}
    limits::HardwareLimits
    gradient::G
    receiver::RX
    transmitter::TX
end

function Scanner(;
    limits=HardwareLimits(),
    gradient=LinearXYZ(),
    receiver=UniformCoilSens(),
    transmitter=UniformTransmit(),
)
    return Scanner(limits, gradient, receiver, transmitter)
end

get_sens(sys::Scanner, x, y, z) = get_sens(sys.receiver, x, y, z)

# Uniform Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, ::Nothing, ::UniformCoilSens, Mxy)
    sig .= @views transpose(sum(Mxy[:, findall(seq.ADC[2:end])]; dims=1))
    return nothing
end

# Uniform other sim methods
function acquire_signal!(sig, sample, ::UniformCoilSens, Mxy)
    sig[sample, :] .= sum(Mxy)
    return nothing
end

# Birdcage Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, ::Nothing, p, receiver::BirdcageCoilSens, Mxy::AbstractMatrix)
    sig .= transpose(Mxy) * get_sens(receiver, p.x, p.y, p.z)
    return nothing
end

# Birdcage other sim methods with motion
function acquire_signal!(sig, sample, p, receiver::BirdcageCoilSens, Mxy::AbstractMatrix)
    sig[sample, :] .= @views vec(transpose(Mxy[:, sample]) * get_sens(receiver, p.x[:, sample], p.y[:, sample], p.z[:, sample]))
    return nothing
end

# Birdcage other sim methods with no motion
function acquire_signal!(sig, sample, p, receiver::BirdcageCoilSens, Mxy::AbstractVector)
    sig[sample, :] .= vec(transpose(Mxy) * get_sens(receiver, p.x, p.y, p.z))
    return nothing
end

# Arbitary Bloch Simple run_spin_precession not sample based
function acquire_signal!(sig, ::Nothing, p, receiver::ArbitraryRFRxCoils, Mxy::AbstractMatrix)
    interpolated_coil_sens = similar(receiver.coil_sens, length(p.x), size(receiver.coil_sens, 4))
    for i in 1:size(receiver.coil_sens, 4)
        base_itp = GriddedInterpolation((receiver.x, receiver.y, receiver.z), receiver.coil_sens[:,:,:,i], Gridded(Linear()))
        itp = extrapolate(base_itp, 0f0)
        interpolated_coil_sens[:,i] = itp.(p.x, p.y, p.z)
        sig[:, i] .= vec(sum(interpolated_coil_sens[:, i] .* Mxy, dims=1))
    end
    return nothing
end

# Arbitary other sim methods with motion
function acquire_signal!(sig, sample, p, receiver::ArbitraryRFRxCoils, Mxy::AbstractMatrix)
    interpolated_coil_sens = similar(receiver.coil_sens, length(p.x), size(receiver.coil_sens, 4))
    for i in 1:size(receiver.coil_sens, 4)
        base_itp = GriddedInterpolation((receiver.x, receiver.y, receiver.z), receiver.coil_sens[:,:,:,i], Gridded(Linear()))
        itp = extrapolate(base_itp, 0f0)
        interpolated_coil_sens[:,i] = itp.(p.x[:, sample], p.y[:, sample], p.z[:, sample])
        sig[sample, i] = sum(interpolated_coil_sens[:, i] .* Mxy)
    end
    return nothing
end

# Arbitary other sim methods with no motion
function acquire_signal!(sig, sample, p, receiver::ArbitraryRFRxCoils, Mxy::AbstractVector)
    interpolated_coil_sens = similar(receiver.coil_sens, length(p.x), size(receiver.coil_sens, 4))
    for i in 1:size(receiver.coil_sens, 4)
        base_itp = GriddedInterpolation((receiver.x, receiver.y, receiver.z), receiver.coil_sens[:,:,:,i], Gridded(Linear()))
        itp = extrapolate(base_itp, 0f0)
        interpolated_coil_sens[:,i] = itp.(p.x, p.y, p.z)
        sig[sample, i] = sum(interpolated_coil_sens[:, i] .* Mxy)
    end
    return nothing
end

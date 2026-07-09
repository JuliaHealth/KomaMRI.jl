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
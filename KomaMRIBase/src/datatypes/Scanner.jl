include("scanner/GradientCoils.jl")
include("scanner/RFTransmitSystem.jl")
include("scanner/RFReceiveSystem.jl")
include("scanner/HardwareLimits.jl")

export HardwareLimits, Scanner,
    AbstractGradientSystem, LinearXYZ,
    AbstractRFReceiveSystem, UniformCoilSens, BirdcageCoilSens, ArbitraryCoilSens,
    AbstractRFTransmitSystem, UniformTransmit,
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
Base.@kwdef struct Scanner{
    G<:AbstractGradientSystem,
    RX<:AbstractRFReceiveSystem,
    TX<:AbstractRFTransmitSystem,
}
    limits::HardwareLimits = HardwareLimits()
    gradient::G = LinearXYZ()
    receiver::RX = UniformCoilSens()
    transmitter::TX = UniformTransmit()
end

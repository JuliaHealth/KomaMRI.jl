"""
Supertype for receive-coil sensitivity models.

Subtypes must implement [`get_n_coils`](@ref) and [`get_sens`](@ref).
"""
abstract type AbstractRFReceiveSystem end

"""Return the number of receive channels in `receiver`."""
function get_n_coils end

"""
    get_sens(receiver, x, y, z)

Evaluate the complex receive sensitivities at the supplied spin positions.
The result has size `(length(x), get_n_coils(receiver))`.
"""
function get_sens end
function get_sens! end

function get_sens!(sens, receiver::AbstractRFReceiveSystem, x, y, z)
    copyto!(sens, get_sens(receiver, x, y, z))
    return sens
end

include("UniformCoilSens.jl")
include("BirdcageCoilSens.jl")
include("ArbitraryCoilSens.jl")

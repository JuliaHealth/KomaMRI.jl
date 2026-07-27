"""Supertype for RF transmit-system models."""
abstract type AbstractRFTransmitSystem end

"""Spatially uniform RF transmit model."""
struct UniformTransmit <: AbstractRFTransmitSystem end

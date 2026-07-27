"""Supertype for scanner gradient-system models."""
abstract type AbstractGradientSystem end

"""Ideal linear gradients along the Cartesian x, y, and z axes."""
struct LinearXYZ <: AbstractGradientSystem end

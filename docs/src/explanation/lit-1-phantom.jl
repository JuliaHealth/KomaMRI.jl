# # Phantom

using KomaMRI # hide

# The first input argument that **KomaMRI** needs for simulating is the phantom. 

# This section goes over the concept of digital phantom 
# and shows how it applies to the specific case of **KomaMRI**. 
# We'll go into detail about the [`Phantom`](@ref) structure 
# and present the ''.phantom'' file format, which makes it easy 
# to share phantoms and reproduce experiments on any computer.

# ## Digital Phantom

# A digital phantom is basically a computer model of a physical object 
# (like the human body or a body part) which is used in simulations 
# to mimic the characteristics and behaviour that would be obtained 
# from real MRI. Instead of using a physical object for testing, 
# the digital phantom allows for virtual experiments.

# This computer model should essentially contain information about 
# the position and/or displacements of the tissues, as well as 
# their MRI-related (T1, T2, PD, off-resonance...) values.

# ## KomaMRI Phantom Overview

# **KomaMRI** relies on the [`Phantom`](@ref) struct to define its digital phantom:
# ```julia
# @with_kw mutable struct Phantom{T<:Real}
#     name::String = "spins"
#     x::AbstractVector{T}
#     y::AbstractVector{T}   = zeros(eltype(x), size(x))
#     z::AbstractVector{T}   = zeros(eltype(x), size(x))
#     ρ::AbstractVector{T}   = ones(eltype(x), size(x))
#     T1::AbstractVector{T}  = ones(eltype(x), size(x)) * 1_000_000
#     T2::AbstractVector{T}  = ones(eltype(x), size(x)) * 1_000_000
#     T2s::AbstractVector{T} = ones(eltype(x), size(x)) * 1_000_000
#     #Off-resonance related
#     Δw::AbstractVector{T}  = zeros(eltype(x), size(x))
#     #Diffusion
#     Dλ1::AbstractVector{T} = zeros(eltype(x), size(x))
#     Dλ2::AbstractVector{T} = zeros(eltype(x), size(x))
#     Dθ::AbstractVector{T}  = zeros(eltype(x), size(x))
#     #Motion
#     motion::AbstractMotion{T} = NoMotion{eltype(x)}() 
# end
# ```

## Phantom File Format
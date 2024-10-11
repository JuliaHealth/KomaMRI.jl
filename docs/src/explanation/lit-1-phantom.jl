# # Phantom

using KomaMRI #hide

# The first input argument that **KomaMRI** needs for simulating is the phantom. 

# This section goes over the concept of digital phantom and shows how it applies to the specific case
# of **KomaMRI**. We'll go into detail about the [`Phantom`](@ref) structure and its supported operations.

# ## Digital Phantom

# A digital phantom is basically a computer model of a physical object (like the human body or a body part) 
# which is used  in simulations to mimic the characteristics and behaviour that would be obtained from
# real MRI. Instead of using a physical object for testing, the digital phantom allows for virtual experiments.

# This computer model should essentially contain information about the position and/or displacements
# of the tissues, as well as their MRI-related (T1, T2, PD, off-resonance...) values. 

# ## KomaMRI Phantom Overview
# In Koma, a phantom is made up of a set of spins (which in many cases are also known as ''isochromats'').
# Each spin is independent of the others in terms of properties, position and state.
# This is a key feature of **KomaMRI**, as it is explained in the [Simulation](6-simulation.md) section.

# Let's take a look at the definition of the [`Phantom`](@ref) struct 
# inside Koma's source code to see what it looks like:
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

# This structure consists of several elements. Most of them are vectors, except for
# the `name` (self-explanatory) and `motion` (explained below) fields.
# These vectors represent object properties, with each element holding a value associated 
# with a single magnetization (i.e. a single spin).
# Specifically, `x`, `y` and `z` are the spatial (starting) coordinates of each spin. 
# `ρ` stands for the proton density, and `T1`, `T2` and `T2s` (standing for T2*) 
# are the well-known relaxation times. `Δw` accounts for off-resonance effects.
# `Dλ1`, `Dλ2` and `Dθ` are diffusion-related fields which are not in use at the moment.
# Last, the `motion` field stands for spin displacements, which are added to `x`, `y` and `z`
# when simulating in order to obtain the spin positions at each time step. For more information about
# motion, refer to [Motion](2-motion.md) section.

# To get an even better understanding on how it works, let's look at an example of a brain phantom:

obj = brain_phantom2D()
# ```julia-repl
# Phantom{Float64}
#   name: String "brain2D_axial"
#   x: Array{Float64}((6506,)) [-0.084, -0.084, -0.084, -0.084, -0.084, -0.084, -0.084, -0.084, -0.084, -0.084  …  0.084, 0.084, 0.084, 0.084, 0.086, 0.086, 0.086, 0.086, 0.086, 0.086]
#   y: Array{Float64}((6506,)) [-0.03, -0.028, -0.026, -0.024, -0.022, -0.02, -0.018, -0.016, -0.014, -0.012  …  0.006, 0.008, 0.01, 0.012, -0.008, -0.006, -0.004, -0.002, 0.0, 0.002]
#   z: Array{Float64}((6506,)) [-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0  …  0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0]
#   ρ: Array{Float64}((6506,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  …  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#   T1: Array{Float64}((6506,)) [0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569  …  0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569, 0.569]
#   T2: Array{Float64}((6506,)) [0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329  …  0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329, 0.329]
#   T2s: Array{Float64}((6506,)) [0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058  …  0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058, 0.058]
#   Δw: Array{Float64}((6506,)) [-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0  …  -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0]
#   Dλ1: Array{Float64}((6506,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#   Dλ2: Array{Float64}((6506,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#   Dθ: Array{Float64}((6506,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#   motion: NoMotion{Float64} NoMotion{Float64}()
# ```

# You can visualize the **Phantom** struct using the [`plot_phantom_map`](@ref) function, 
# which is part of the **KomaMRIPlots** subdependency. This function plots the magnitude of a property for 
# each magnetization at a specific spatial position. You can observe properties such as proton density 
# and relaxation times, so feel free to replace the `:T1` symbol with another property of the phantom in the example below:

p1 = plot_phantom_map(obj, :T1; height=450)

#md savefig(p1, "../assets/doc-1-phantom.html") #hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/doc-1-phantom.html" style="width:85%; height:470px;"></object></center>
#md # ```

# You can access and filter information for the all the field names of a **Phantom** using the dot notation:

# ```julia-repl
# julia> obj.name
# "brain2D_axial"
# ```

# ```julia-repl
# julia> obj.x
# 6506-element Vector{Float64}:
#  -0.084
#  -0.084
#  -0.084
#   ⋮
#   0.086
#   0.086
#   0.086
# ```

# ```julia-repl
# julia> obj.motion
# NoMotion{Float64}()
# ```

# ## Phantom Operations

# In addition, **KomaMRI** supports some phantom operations:

# ### Phantom Subset

# It is possible to access a subset of spins in a **Phantom** by slicing or indexing. The result will also be a 
# **Phantom** struct, allowing you to perform the same operations as you would with a full Phantom:

obj[1:2000]
p2 = plot_phantom_map(obj[1:1000], :T2 ; height=450) #hide

#md savefig(p2, "../assets/tut-5-phantom-subset.html") #hide
#jl display(p2)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/tut-5-phantom-subset.html" style="width:85%; height:470px;"></object></center>
#md # ```

# ### Combination of Phantoms

# In the same way, we can add two or more phantoms, resulting in another [`Phantom`](@ref) struct:
obj2 = pelvis_phantom2D()
obj2.x .+= 0.1; obj.x.-= 0.1 #hide
obj_sum = obj + obj2
p3 = plot_phantom_map(obj_sum, :T1 ; height=450) #hide

#md savefig(p3, "../assets/tut-5-phantom-sum.html") #hide
#jl display(p3)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/tut-5-phantom-sum.html" style="width:85%; height:470px;"></object></center>
#md # ```

# ### Scalar multiplication of a Phantom

# Finally, multiplying a phantom by a scalar multiplies its proton density (`ρ`) by that amount:
obj_mul = 3*obj

# ```julia-repl
# julia> obj.ρ
# 6506-element Vector{Float64}:
#  1.0
#  1.0
#  1.0
#  ⋮
#  1.0
#  1.0
#  1.0
#  
# julia> obj_mul.ρ
# 6506-element Vector{Float64}:
#  3.0
#  3.0
#  3.0
#  ⋮
#  3.0
#  3.0
#  3.0
# ```

# ## Phantom Storage and Sharing
# Phantoms can be stored and shared thanks to our new [Phantom File Format](3-phantom-format.md).
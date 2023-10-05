# # Gradient Design

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./FILE_NAME.jl)
#md # [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./FILE_NAME.ipynb)

# Before doing anything, let's import the KomaMRI package

using KomaMRI

# ## Gradient Trapezoidal Waveform

# The most simple definition of a gradient is a trapezoidal waveform. Let's create a sequence
# (a sequence with just one block), with 3 gradients:

## Define parameters and gradients
## A: amplitude of the gradient
## T: duration of the flat top
## ζ: rise and fall duration
A, T, ζ = 10e-3, 1e-3, 0.2e-3
gx = Grad(   A, T, ζ)
gy = Grad( 2*A, T, ζ)
gz = Grad(.5*A, T, ζ)

## Define the sequence (just one block) with the 3 gradients
grs = reshape([gx; gy; gz], :, 1)
seq = Sequence(grs)

# To visually check that everithing looks ok, let's plot the sequence with the defined
# trapezoidal gradients

p1 = plot_seq(seq; slider=false, height=300)
#md savefig(p1, "../../assets/FOLDER_NAME/1-gr-trap.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/1-gr-trap.html" style="width:100%; height:320px;"></object>
#md # ```


# ## Gradient Uniformly-Sampled Waveform

# You can also create gradients waveforms with more flexible shapes. Let's create a gradient
# in one axis, with an arbritary amplitude shape at evenly spaced times:

## Define parameters
## T: duration of the arbitrary waveform
## ζ: rise and fall duration
## A: vector with the amplitudes of the gradient
T, ζ = 1e-3, 0.05e-3
A = 10e-3 * (sin.(0.:2π/10:2π) .+ 1)

## Define the gradient and the sequence (with just one block)
gx = Grad(A, T, ζ)
seq = Sequence([gx])

p2 = plot_seq(seq; slider=false, height=300)
#md savefig(p2, "../../assets/FOLDER_NAME/2-gr-wave.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/2-gr-wave.html" style="width:100%; height:320px;"></object>
#md # ```

# Note that the waveform always defines rise and fall times, so the gradients always starts
# and ends at zero amplitude.


# ## Gradient Time-Shaped Waveform

# You can also create a waveform sampled at arbitraty time positions (as long as they are in
# increasing order). Let's create such gradient:

## Define parameters
## A: vector with the amplitudes of the gradient
## T: vector with the distances between two consecutive time samples
## ζ: rise and fall duration

A = 10e-3 * (sin.(0.:π/10:π) .+ 1)
T = 10e-3 * [((-1:0.2:-0.2).^2); ((0.2:0.2:1).^2)]
ζ = 0.05e-3

## Define the gradient and the sequence (with just one block)
gx = Grad(A, T, ζ)
seq = Sequence([gx])

p3 = plot_seq(seq; slider=false, height=300)
#md savefig(p3, "../../assets/FOLDER_NAME/3-gr-wave.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/3-gr-wave.html" style="width:100%; height:320px;"></object>
#md # ```

# Note that in this case the length of the amplitud vector must have one sample more than the
# samples of the vector of time distances.

# # GR Design

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./FILE_NAME.jl)
#md # [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./FILE_NAME.ipynb)

# Let's import the KomaMRI package

using KomaMRI

# ## GR Trapezoidal definition

# The most simple definition of a gradient is a trapezoidal waveform. Let's create a sequence
# (a sequence with just one block), with 3 x,y,z gradients:

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

seq.GR[1,1]
seq.GR[1,1].A
seq.GR[1,1].T
seq.GR[1,1].rise
seq.GR[1,1].fall

seq.GR[2,1]
seq.GR[2,1].A
seq.GR[2,1].T
seq.GR[2,1].rise
seq.GR[2,1].fall

seq.GR[3,1]
seq.GR[3,1].A
seq.GR[3,1].T
seq.GR[3,1].rise
seq.GR[3,1].fall

# ## GR Arbitrary waveform

A = 10e-3 * [1; .65; .2; .1; .2; .65; 1]
T, ζ = 1e-3, 0.2e-3
gx = Grad(   A, T, ζ)
gy = Grad( 2*A, T, ζ)
gz = Grad(.5*A, T, ζ)
gr_m = reshape([gx; gy; gz], :, 1)
seq = Sequence(gr_m)

p2 = plot_seq(seq; slider=false, height=300)
#md savefig(p2, "../../assets/FOLDER_NAME/2-gr-wave.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/2-gr-wave.html" style="width:100%; height:320px;"></object>
#md # ```

# ## GR x y z concatenation

# ## GR time concatenation

# ## GR rotation

# ## GR example: moment compensated diffusion

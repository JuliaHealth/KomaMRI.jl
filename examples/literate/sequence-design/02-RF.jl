# # RF Design

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./FILE_NAME.jl)
#md # [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./FILE_NAME.ipynb)

# Let's import the KomaMRI package

using KomaMRI

# ## RF hard definition

# The first block containing an RF pulse with a flip-angle of 90 deg,
# To define an RF object, at least we need to add define its amplitude and duration.
# In paticular for a RF-hard definition (flip-angle of 90 deg), the duration of the RF object
# is a function of it's amplitude

A = 2e-6                    # 2 uT RF amplitude
T = π / 2 / (2π * γ * A)    # required duration for a 90 deg RF pulse
rf = RF(A, T)

# The previous definition of the RF would be enough, however in order to work with it, we
# need to add this RF object into a Sequence object. The minimal definition for a Sequence
# with an RF requires Gradient and RF matrices (3x1 and 1x1 dimensions respectively),
# thus the Sequence with an RF-hard pulse can be defined like so:

g0 = Grad(0, 0)
gr_m = reshape([g0; g0; g0], :, 1)
rf_m = reshape([rf], :, 1)
seq = Sequence(gr_m, rf_m)

# We can plot the definition of the RF-hard signal by ploting the single block sequence:

p1 = plot_seq(seq; slider=false, height=300)
#md savefig(p1, "../../assets/FOLDER_NAME/1-rf-hard.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/1-rf-hard.html" style="width:100%; height:320px;"></object>
#md # ```

# Note that the RF object can be extracted by manipulating the "first RF" of the sequence
# like so:

seq.RF[1]
seq.RF[1].A
seq.RF[1].T

# ## RF waveform

# We can also define arbitrary waveforms for the RF object, this can be done considering
# a vector samples for the RF amplitude and asuming equally spaced samples for a
# duration of time:

A = 10e-6 * [.25; -.5; 1.; -.5; .25]
T = 5.
rf = RF(A, T)

# Then work with a single block sequence:

rf_m = reshape([rf], :, 1)
seq = Sequence(gr_m, rf_m)

# Then visualize de the RF object with waveform created:

p1 = plot_seq(seq; slider=false, height=300)
#md savefig(p1, "../../assets/FOLDER_NAME/2-rf-soft.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/2-rf-soft.html" style="width:100%; height:320px;"></object>
#md # ```

# And fanally manipulate the RF first RF of the sequence:

seq.RF[1]
seq.RF[1].A
seq.RF[1].T

# ## RF time concatenation

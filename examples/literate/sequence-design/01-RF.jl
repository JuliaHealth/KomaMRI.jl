# # RF Design

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./FILE_NAME.jl)
#md # [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./FILE_NAME.ipynb)

# First of all, let's use the KomaMRI package and define the default scanner.

using KomaMRI
sys = Scanner() # default hardware definition

# The free induction decay is the simplest observable NMR signal.
# This signal is the one that follows a single tipping RF pulse.
# To recreate this experiment, we will need to define a `Sequence` with 2 blocks.

# The first block containing an RF pulse with a flip-angle of 90 deg,

ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF,durRF)

# and the second block containing the ADC.

nADC = 8192         # number of acquisition samples
durADC = 250e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay)

# Finally, we concatenate the sequence blocks to create
# the final sequence.

seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
p1 = plot_seq(seq; slider=false, height=300)
#md savefig(p1, "../../assets/FOLDER_NAME/1-seq.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <object type="text/html" data="../../../assets/FOLDER_NAME/1-seq.html" style="width:100%; height:320px;"></object>
#md # ```

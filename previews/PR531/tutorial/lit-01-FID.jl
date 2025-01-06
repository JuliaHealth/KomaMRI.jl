# # Free Induction Decay

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
#md savefig(p1, "../assets/1-seq.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/1-seq.html" style="width:100%; height:320px;"></object>
#md # ```

# Now, we will define a `Phantom` with a single spin at ``x=0``
# with ``T_1=1000\,\mathrm{ms}`` and ``T_2=100\,\mathrm{ms}``.

obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])

# Finally, to simulate we will need to use the function [`simulate`](@ref).

raw = simulate(obj, seq, sys)

# To plot the results we will need to use the [`plot_signal`](@ref) function

p2 = plot_signal(raw; slider=false, height=300)
#md savefig(p2, "../assets/1-signal.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/1-signal.html" style="width:100%; height:320px;"></object>
#md # ```

# Nice!, we can see that ``S(t)`` follows an
# exponential decay ``\exp(-t/T_2)`` as expected.

# For a little bit of spiciness, let's add **off-resonance** to our example.
# We will use ``\Delta f=-100\,\mathrm{Hz}``.
# For this, we will need to add a definition for `Δw` in our `Phantom`

obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3], Δw=[-2π*100])# and simulate again.

raw = simulate(obj, seq, sys)
p3 = plot_signal(raw; slider=false, height=300)
#md savefig(p3, "../assets/1-signal2.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/1-signal2.html" style="width:100%; height:320px;"></object>
#md # ```

# The signal now follows an exponential of the
# form ``\exp(-t/T_2)\cdot\exp(-i\Delta\omega t)``.
# The addition of ``\exp(-i\Delta\omega t)`` to the signal
# will generate a shift in the image space (Fourier shifting property).
# This effect will be better visualized and explained in later examples.

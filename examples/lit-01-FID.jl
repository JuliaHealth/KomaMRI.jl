# # Free Induction Decay

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](@__REPO_ROOT_URL__/examples/lit-01-FID.jl)

filename = last(splitpath(@__FILE__)) # hide
isFileMD = occursin(".md", filename) # hide
isFileJL = occursin(".jl", filename) # hide

#md # First of all, let's use the KomaMRI package and define the default scanner.

using KomaMRI
sys = Scanner() # default hardware definition
nothing # hide

#md # The free induction decay is the simplest observable NMR signal.
#md # This signal is the one that follows a single tipping RF pulse.
#md # To recreate this experiment, we will need to define a `Sequence` with 2 blocks.

#md # The first block containing an RF pulse with a flip-angle of 90 deg,

ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF,durRF)
nothing # hide

#md # and the second block containing the ADC.

nADC = 8192         # number of acquisition samples
durADC = 250e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay)
nothing # hide

#md # Finally, we concatenate the sequence blocks to create
#md # the final sequence (for more info. refer to
#md # [Sequence Structure](useful-information.md#Sequence-Structure)).

seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
p1 = plot_seq(seq; slider=false, height=300)
if isFileMD savefig(p1, "../assets/1-seq.html") end # hide
if isFileJL display(p1) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/1-seq.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Now, we will define a `Phantom` with a single spin at ``x=0``
#md # with ``T_1=1000\,\mathrm{ms}`` and ``T_2=100\,\mathrm{ms}``.

obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])
nothing # hide

#md # Finally, to simulate we will need to use the function [`simulate`](@ref).

raw = simulate(obj, seq, sys)

#md # To plot the results we will need to use the [`plot_signal`](@ref) function

p2 = plot_signal(raw; slider=false, height=300)
if isFileMD savefig(p2, "../assets/1-signal.html") end # hide
if isFileJL display(p2) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/1-signal.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Nice!, we can see that ``S(t)`` follows an
#md # exponential decay ``\exp(-t/T_2)`` as expected.

#md # For a little bit of spiciness, let's add **off-resonance** to our example.
#md # We will use ``\Delta f=-100\,\mathrm{Hz}``.
#md # For this, we will need to add a definition for `Δw` in our `Phantom`

obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3], Δw=[-2π*100])
nothing # hide

#md # and simulate again.

raw = simulate(obj, seq, sys)
p3 = plot_signal(raw; slider=false, height=300)
if isFileMD savefig(p3, "../assets/1-signal2.html") end # hide
if isFileJL display(p3) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/1-signal2.html" style="width:100%; height:320px;"></object>
#md # ```

#md # The signal now follows an exponential of the
#md # form ``\exp(-t/T_2)\cdot\exp(-i\Delta\omega t)``.
#md # The addition of ``\exp(-i\Delta\omega t)`` to the signal
#md # will generate a shift in the image space (Fourier shifting property).
#md # This effect will be better visualized and explained in later examples.

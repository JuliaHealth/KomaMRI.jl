```@meta
EditURL = "<unknown>/../examples/lit/examples/01-FID.jl"
```

# Free Induction Decay

First of all, let's use the KomaMRI package and define the default scanner.

````@example 01-FID
using KomaMRI
sys = Scanner() # default hardware definition
nothing # hide
````

The free induction decay is the simplest observable NMR signal.
This signal is the one that follows a single tipping RF pulse.
To recreate this experiment, we will need to define a `Sequence` with 2 blocks.

The first block containing an RF pulse with a flip-angle of 90 deg,

````@example 01-FID
ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
exc = RF(ampRF,durRF)
nothing # hide
````

and the second block containing the ADC.

````@example 01-FID
nADC = 8192         # number of acquisition samples
durADC = 250e-3     # duration of the acquisition
delay =  1e-3       # small delay
acq = ADC(nADC, durADC, delay)
nothing # hide
````

Finally, we concatenate the sequence blocks to create
the final sequence (for more info. refer to
[Sequence Structure](useful-information.md#Sequence-Structure)).

````@example 01-FID
seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
p = plot_seq(seq; slider=false, height=300)
savefig(p, "../assets/1-seq.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/1-seq.html" style="width:100%; height:320px;"></object>
```

Now, we will define a `Phantom` with a single spin at ``x=0``
with ``T_1=1000\,\mathrm{ms}`` and ``T_2=100\,\mathrm{ms}``.

````@example 01-FID
obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3])
nothing # hide
````

Finally, to simulate we will need to use the function [`simulate`](@ref).

````@example 01-FID
raw = simulate(obj, seq, sys)
````

To plot the results we will need to use the [`plot_signal`](@ref) function

````@example 01-FID
p = plot_signal(raw; slider=false, height=300)
savefig(p, "../assets/1-signal.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/1-signal.html" style="width:100%; height:320px;"></object>
```

Nice!, we can see that ``S(t)`` follows an
exponential decay ``\exp(-t/T_2)`` as expected.

For a little bit of spiciness, let's add **off-resonance** to our example.
We will use ``\Delta f=-100\,\mathrm{Hz}``.
For this, we will need to add a definition for `Δw` in our `Phantom`

````@example 01-FID
obj = Phantom{Float64}(x=[0.], T1=[1000e-3], T2=[100e-3], Δw=[-2π*100])
nothing # hide
````

and simulate again.

````@example 01-FID
raw = simulate(obj, seq, sys)
p = plot_signal(raw; slider=false, height=300)
savefig(p, "../assets/1-signal2.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/1-signal2.html" style="width:100%; height:320px;"></object>
```

The signal now follows an exponential of the
form ``\exp(-t/T_2)\cdot\exp(-i\Delta\omega t)``.
The addition of ``\exp(-i\Delta\omega t)`` to the signal
will generate a shift in the image space (Fourier shifting property).
This effect will be better visualized and explained in later examples.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


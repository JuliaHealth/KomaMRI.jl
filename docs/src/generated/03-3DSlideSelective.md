```@meta
EditURL = "<unknown>/../examples/lit/examples/03-3DSlideSelective.jl"
```

# Slice Selective Acquisition of 3D Phantom

First of all, let's use the KomaMRI package and define the default scanner.

````@example 03-3DSlideSelective
using KomaMRI
sys = Scanner() # default hardware definition
nothing # hide
````

While in the previous examples we simulated using hard RF pulses,
in this demonstration we will illustrate the principles of slice selection.
First, let's import a 3D phantom, in this case a brain slab
(thickness of ``2\,\mathrm{cm}``), by calling the function [`brain_phantom3D`](@ref).

````@example 03-3DSlideSelective
obj = brain_phantom3D()
obj.Δw .= 0 # Removes the off-resonance
p = plot_phantom_map(obj, :T2 ; height=400)
savefig(p, "../assets/3-phantom.html") # hide
nothing # hide
````

```@raw html
<center><object type="text/html" data="../../assets/3-phantom.html" style="width:50%; height:420px;"></object></center>
```

Now, we are going to import a sequence which acquires
3 slices in the longitudinal axis. Note that the sequence
contains three EPIs to acquire 3 slices of the phantom.

```julia
seq = read_seq("examples/1.sequences/epi_multislice.seq")
p = plot_seq(seq; range=[0,10], height=400)
```

````@example 03-3DSlideSelective
seq = read_seq("../../../examples/1.sequences/epi_multislice.seq") # hide
p = plot_seq(seq; range=[0,10], height=400) # hide
savefig(p, "../assets/3-seq.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/3-seq.html" style="width:100%; height:420px;"></object>
```

We can take a look to the slice profiles by using the function [`simulate_slice_profile`](@ref):

````@example 03-3DSlideSelective
z = range(-2., 2., 200) * 1e-2; # -2 to 2 cm
rf1, rf2, rf3 = findall(KomaMRI.is_RF_on.(seq))
M1 = simulate_slice_profile(seq[rf1]; z)
M2 = simulate_slice_profile(seq[rf2]; z)
M3 = simulate_slice_profile(seq[rf3]; z)

using PlotlyJS # hide
p1 = scatter(x=z*1e2, y=abs.(M1.xy), name="Slice 1") # hide
p2 = scatter(x=z*1e2, y=abs.(M2.xy), name="Slice 2") # hide
p3 = scatter(x=z*1e2, y=abs.(M3.xy), name="Slice 3") # hide
p = plot([p1,p2,p3], Layout(xaxis=attr(title="z [cm]"), height=300,margin=attr(t=40,l=0,r=0), title="Slice profiles for the slice-selective sequence")) # hide
savefig(p, "../assets/3-profile.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/3-profile.html" style="width:100%; height:320px;"></object>
```

Now let's simulate the acquisition.
Notice the three echoes, one for every slice excitation.

````@example 03-3DSlideSelective
raw = simulate(obj, seq, sys; simParams=Dict{String,Any}("Nblocks"=>20))
p = plot_signal(raw; slider=false, height=300)
savefig(p, "../assets/3-signal.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/3-signal.html" style="width:100%; height:320px;"></object>
```

Finally, we reconstruct the acquiered images.

````@example 03-3DSlideSelective
# Get the acquisition data
acq = AcquisitionData(raw)

# Setting up the reconstruction parameters and perform reconstruction
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

# Plotting the slices
p1 = plot_image(abs.(image[:, :, 1]); height=360, title="Slice 1")
p2 = plot_image(abs.(image[:, :, 2]); height=360, title="Slice 2")
p3 = plot_image(abs.(image[:, :, 3]); height=360, title="Slice 3")
savefig(p1, "../assets/3-recon1.html") # hide
savefig(p2, "../assets/3-recon2.html") # hide
savefig(p3, "../assets/3-recon3.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../assets/3-recon1.html" style="width:50%; height:380px;"></object><object type="text/html" data="../../assets/3-recon2.html" style="width:50%; height:380px;"></object>
```
```@raw html
<center><object type="text/html" data="../../assets/3-recon3.html" style="width:50%; height:380px;"></object></center>
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


```@meta
EditURL = "<unknown>/lit/examples/example-02.jl"
```

# Chemical Shift in an EPI sequence

First of all, let's use the KomaMRI package and define the default scanner.

````@example example-02
using KomaMRI
sys = Scanner() # default hardware definition
nothing # hide
````

For a more realistic example, we will use a brain phantom.

````@example example-02
obj = brain_phantom2D() # a slice of a brain
p1 = plot_phantom_map(obj, :T2 ; height=400)
p2 = plot_phantom_map(obj, :Δw ; height=400)
savefig(p1, "../../assets/2-phantom1.html") # hide
savefig(p2, "../../assets/2-phantom2.html") # hide
nothing # hide
````

At the left, you can see the ``T_2`` map of the phantom,
and at the right, the off-resonance ``\Delta\omega``.
In this example, the fat is the only source of off-resonance
(with ``\Delta f =  -220\,\mathrm{Hz}``) and you can see
it in black in the off-resonance map.

```@raw html
<object type="text/html" data="../../../assets/2-phantom1.html" style="width:50%; height:420px;"></object><object type="text/html" data="../../../assets/2-phantom2.html" style="width:50%; height:420px;"></object>
```

Then, we will load an EPI sequence, that is well known
for being affected by off-resonance. With this sequence,
we will be able visualize the effect of the chemical shift.

```julia
seq = read_seq("examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq")
p = plot_seq(seq; range=[0 40], slider=true, height=300)
```
```@setup example-02
seq = read_seq("../../../../examples/3.koma_paper/comparison_jemris/sequences/EPI/epi_100x100_TE100_FOV230.seq")
```

````@example example-02
savefig(plot_seq(seq; range=[0 40], slider=true, height=300), "../../assets/2-seq.html") # hide
nothing # hide
````

Feel free to explore the sequence's plot 🔍 below!

```@raw html
<object type="text/html" data="../../../assets/2-seq.html" style="width:100%; height:320px;"></object>
```

If we simulate this sequence we will end up with the following signal.

```julia
raw = simulate(obj, seq, sys)
p = plot_signal(raw; range=[98.4 103.4] , height=300)
```
```@setup example-02
raw = simulate(obj, seq, sys)
```

````@example example-02
savefig(plot_signal(raw; range=[98.4 103.4] , height=300), "../../assets/2-signal.html") # hide
nothing # hide
````

```@raw html
<object type="text/html" data="../../../assets/2-signal.html" style="width:100%; height:320px;"></object>
```

Now, we need to inspect what effect the off-resonance
had in the reconstructed image. As you can see,
the fat layer is now shifted to a different position 🤯,
this is why the effect is called chemical shift!

````@example example-02
# Get the acquisition data
acq = AcquisitionData(raw)
acq.traj[1].circular = false #This is to remove a circular mask

# Setting up the reconstruction parameters
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

# Plotting the recon
slice_abs = abs.(image[:, :, 1])
p = plot_image(slice_abs; height=400)
savefig(p, "../../assets/2-recon.html") # hide
nothing # hide
````

```@raw html
<center><object type="text/html" data="../../../assets/2-recon.html" style="width:65%; height:420px;"></object></center>
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


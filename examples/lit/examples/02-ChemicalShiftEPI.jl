# # Chemical Shift in an EPI sequence

#md # First of all, let's use the KomaMRI package and define the default scanner.

using KomaMRI
sys = Scanner() # default hardware definition
nothing # hide

#md # For a more realistic example, we will use a brain phantom.

obj = brain_phantom2D() # a slice of a brain
p1 = plot_phantom_map(obj, :T2 ; height=400)
p2 = plot_phantom_map(obj, :Δw ; height=400)
savefig(p1, "../assets/2-phantom1.html") # hide
savefig(p2, "../assets/2-phantom2.html") # hide
nothing # hide

#md # At the left, you can see the ``T_2`` map of the phantom,
#md # and at the right, the off-resonance ``\Delta\omega``.
#md # In this example, the fat is the only source of off-resonance
#md # (with ``\Delta f =  -220\,\mathrm{Hz}``) and you can see
#md # it in black in the off-resonance map.

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-phantom1.html" style="width:50%; height:420px;"></object><object type="text/html" data="../../assets/2-phantom2.html" style="width:50%; height:420px;"></object>
#md # ```

#md # Then, we will load an EPI sequence, that is well known
#md # for being affected by off-resonance. With this sequence,
#md # we will be able visualize the effect of the chemical shift.

# ```julia
# seq = read_seq("examples/3.koma_paper/comparison/sequences/EPI/epi_100x100_TE100_FOV230.seq")
# p = plot_seq(seq; range=[0 40], slider=true, height=300)
# ```

seq = read_seq("../../../examples/3.koma_paper/comparison_jemris/sequences/EPI/epi_100x100_TE100_FOV230.seq") # hide
p = plot_seq(seq; range=[0 40], slider=true, height=300) # hide
savefig(p, "../assets/2-seq.html") # hide
nothing # hide

#md # Feel free to explore the sequence's plot 🔍 below!

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-seq.html" style="width:100%; height:320px;"></object>
#md # ```

#md # If we simulate this sequence we will end up with the following signal.

raw = simulate(obj, seq, sys)
p = plot_signal(raw; range=[98.4 103.4] , height=300)
savefig(p, "../assets/2-signal.html") # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-signal.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Now, we need to inspect what effect the off-resonance
#md # had in the reconstructed image. As you can see,
#md # the fat layer is now shifted to a different position 🤯,
#md # this is why the effect is called chemical shift!

## Get the acquisition data
acq = AcquisitionData(raw)
acq.traj[1].circular = false #This is to remove a circular mask

## Setting up the reconstruction parameters
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

## Plotting the recon
slice_abs = abs.(image[:, :, 1])
p = plot_image(slice_abs; height=400)
savefig(p, "../assets/2-recon.html") # hide
nothing # hide

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/2-recon.html" style="width:65%; height:420px;"></object></center>
#md # ```

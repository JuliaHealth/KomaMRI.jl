# # Chemical Shift in an EPI sequence

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](@__REPO_ROOT_URL__/examples/lit-03-ChemicalShiftEPI.jl)

filename = last(splitpath(@__FILE__)) # hide
isFileMD = occursin(".md", filename) # hide
isFileJL = occursin(".jl", filename) # hide

using KomaMRI # hide
sys = Scanner() # hide

#md # For a more realistic example, we will use a brain phantom.

obj = brain_phantom2D() # a slice of a brain
p1 = plot_phantom_map(obj, :T2 ; height=400)
p2 = plot_phantom_map(obj, :Δw ; height=400)
if isFileMD savefig(p1, "../assets/2-phantom1.html") end # hide
if isFileMD savefig(p2, "../assets/2-phantom2.html") end # hide
if isFileJL display(p1) end # hide
if isFileJL display(p2) end # hide
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

seqFile = joinpath(dirname(pathof(KomaMRI)), "../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = read_seq(seqFile)
p3 = plot_seq(seq; range=[0 40], slider=true, height=300)
if isFileMD savefig(p3, "../assets/2-seq.html") end # hide
if isFileJL display(p3) end # hide
nothing # hide

#md # Feel free to explore the sequence's plot 🔍 below!

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-seq.html" style="width:100%; height:320px;"></object>
#md # ```

#md # If we simulate this sequence we will end up with the following signal.

raw = simulate(obj, seq, sys)
p4 = plot_signal(raw; range=[98.4 103.4] , height=300)
if isFileMD savefig(p4, "../assets/2-signal.html") end # hide
if isFileJL display(p4) end # hide
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
p5 = plot_image(slice_abs; height=400)
if isFileMD savefig(p5, "../assets/2-recon.html") end # hide
if isFileJL display(p5) end # hide
nothing # hide

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/2-recon.html" style="width:65%; height:420px;"></object></center>
#md # ```

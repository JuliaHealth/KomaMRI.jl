# # Slice Selective Acquisition of 3D Phantom

#md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](@__REPO_ROOT_URL__/examples/lit-04-3DSliceSelective.jl)

filename = basename(@__FILE__) # hide
isFileMD = occursin(".md", filename) # hide
isFileJL = occursin(".jl", filename) # hide

using KomaMRI # hide
sys = Scanner() # hide

#md # While in the previous examples we simulated using hard RF pulses,
#md # in this demonstration we will illustrate the principles of slice selection.
#md # First, let's import a 3D phantom, in this case a brain slab
#md # (thickness of ``2\,\mathrm{cm}``), by calling the function [`brain_phantom3D`](@ref).

obj = brain_phantom3D()
obj.Δw .= 0 # Removes the off-resonance
p1 = plot_phantom_map(obj, :T2 ; height=400)
if isFileMD savefig(p1, "../assets/3-phantom.html") end # hide
if isFileJL display(p1) end # hide
nothing # hide

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/3-phantom.html" style="width:50%; height:420px;"></object></center>
#md # ```

#md # Now, we are going to import a sequence which acquires
#md # 3 slices in the longitudinal axis. Note that the sequence
#md # contains three EPIs to acquire 3 slices of the phantom.

seqFile = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/epi_multislice.seq")
seq = read_seq(seqFile)
p2 = plot_seq(seq; range=[0,10], height=400)
if isFileMD savefig(p2, "../assets/3-seq.html") end # hide
if isFileJL display(p2) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-seq.html" style="width:100%; height:420px;"></object>
#md # ```

#md # We can take a look to the slice profiles by using the function [`simulate_slice_profile`](@ref):

z = range(-2., 2., 200) * 1e-2; # -2 to 2 cm
rf1, rf2, rf3 = findall(KomaMRI.is_RF_on.(seq))
M1 = simulate_slice_profile(seq[rf1]; z)
M2 = simulate_slice_profile(seq[rf2]; z)
M3 = simulate_slice_profile(seq[rf3]; z)

using PlotlyJS # hide
pa = scatter(x=z*1e2, y=abs.(M1.xy), name="Slice 1") # hide
pb = scatter(x=z*1e2, y=abs.(M2.xy), name="Slice 2") # hide
pc = scatter(x=z*1e2, y=abs.(M3.xy), name="Slice 3") # hide
pd = plot([pa,pb,pc], Layout(xaxis=attr(title="z [cm]"), height=300,margin=attr(t=40,l=0,r=0), title="Slice profiles for the slice-selective sequence")) # hide
if isFileMD savefig(pd, "../assets/3-profile.html") end # hide
if isFileJL display(pd) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-profile.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Now let's simulate the acquisition.
#md # Notice the three echoes, one for every slice excitation.

raw = simulate(obj, seq, sys; simParams=Dict{String,Any}("Nblocks"=>20))
p3 = plot_signal(raw; slider=false, height=300)
if isFileMD savefig(p3, "../assets/3-signal.html") end # hide
if isFileJL display(p3) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-signal.html" style="width:100%; height:320px;"></object>
#md # ```

#md # Finally, we reconstruct the acquiered images.

## Get the acquisition data
acq = AcquisitionData(raw)

## Setting up the reconstruction parameters and perform reconstruction
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

## Plotting the slices
p4 = plot_image(abs.(image[:, :, 1]); height=360, title="Slice 1")
p5 = plot_image(abs.(image[:, :, 2]); height=360, title="Slice 2")
p6 = plot_image(abs.(image[:, :, 3]); height=360, title="Slice 3")
if isFileMD savefig(p4, "../assets/3-recon1.html") end # hide
if isFileMD savefig(p5, "../assets/3-recon2.html") end # hide
if isFileMD savefig(p6, "../assets/3-recon3.html") end # hide
if isFileJL display(p4) end # hide
if isFileJL display(p5) end # hide
if isFileJL display(p6) end # hide
nothing # hide

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-recon1.html" style="width:50%; height:380px;"></object><object type="text/html" data="../../assets/3-recon2.html" style="width:50%; height:380px;"></object>
#md # ```
#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/3-recon3.html" style="width:50%; height:380px;"></object></center>
#md # ```

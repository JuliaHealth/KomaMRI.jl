# # Chemical Shift in an EPI sequence

using KomaMRI # hide
sys = Scanner() # hide

# For a more realistic example, we will use a brain phantom.

obj = brain_phantom2D() # a slice of a brain
p1 = plot_phantom_map(obj, :T2 ; height=400, width=400, view_2d=true)
p2 = plot_phantom_map(obj, :Δw ; height=400, width=400, view_2d=true)
#md savefig(p1, "../assets/2-phantom1.html") # hide
#md savefig(p2, "../assets/2-phantom2.html") # hide
#jl display(p1)
#jl display(p2)

# At the left, you can see the ``T_2`` map of the phantom,
# and at the right, the off-resonance ``\Delta\omega``.
# In this example, the fat is the only source of off-resonance
# (with ``\Delta f =  -220\,\mathrm{Hz}``) and you can see
# it in black in the off-resonance map.

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-phantom1.html" style="width:50%; height:420px;"></object><object type="text/html" data="../../assets/2-phantom2.html" style="width:50%; height:420px;"></object>
#md # ```

# Then, we will load an EPI sequence, that is well known
# for being affected by off-resonance. With this sequence,
# we will be able visualize the effect of the chemical shift.

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = read_seq(seq_file)
p3 = plot_seq(seq; range=[0 40], slider=true, height=300)
#md savefig(p3, "../assets/2-seq.html") # hide
#jl display(p3)

#md # Feel free to explore the sequence's plot 🔍 below!

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-seq.html" style="width:100%; height:320px;"></object>
#md # ```

# If we simulate this sequence we will end up with the following signal.

raw = simulate(obj, seq, sys)
p4 = plot_signal(raw; range=[98.4 103.4] , height=300)
#md savefig(p4, "../assets/2-signal.html") # hide
#jl display(p4)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/2-signal.html" style="width:100%; height:320px;"></object>
#md # ```

# Now, we need to inspect what effect the off-resonance
# had in the reconstructed image. As you can see,
# the fat layer is now shifted to a different position 🤯,
# this is why the effect is called chemical shift!

## Get the acquisition data
acq = AcquisitionData(raw)
acq.traj[1].circular = false #This is to remove the circular mask

## Setting up the reconstruction parameters
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

## Plotting the recon
slice_abs = abs.(image[:, :, 1])
p5 = plot_image(slice_abs; height=400)
#md savefig(p5, "../assets/2-recon.html") # hide
#jl display(p5)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/2-recon.html" style="width:65%; height:420px;"></object></center>
#md # ```

# # Chemical Shift in an EPI sequence

using KomaMRI, Suppressor #hide
sys = Scanner(); #hide

# For a more realistic example, we will use a brain phantom.

obj = brain_phantom2D() # a slice of a brain
p1 = plot_phantom_map(obj, :T2 ; height=400, width=400, view_2d=true)
p2 = plot_phantom_map(obj, :Δw ; height=400, width=400, view_2d=true)
#md [p1 p2] #hide
#jl display([p1 p2])

# At the left, you can see the ``T_2`` map of the phantom,
# and at the right, the off-resonance ``\Delta\omega``.
# In this example, the fat is the only source of off-resonance
# (with ``\Delta f =  -220\,\mathrm{Hz}``) and you can see
# it in black in the off-resonance map.

# Then, we will load an EPI sequence, that is well known
# for being affected by off-resonance. With this sequence,
# we will be able visualize the effect of the chemical shift.

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = @suppress read_seq(seq_file)
p3 = plot_seq(seq; range=[0 40], slider=true, height=300)
#jl display(p3);

#md # Feel free to explore the sequence's plot 🔍 below!

# If we simulate this sequence we will end up with the following signal.

raw = @suppress simulate(obj, seq, sys)
p4 = plot_signal(raw; range=[98.4 103.4] , height=300)
#jl display(p4);

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
#jl display(p5);

# # Slice-Selective Acquisition of 3D Phantom

using KomaMRI # hide
sys = Scanner() # hide

# While in the previous examples we simulated using hard RF pulses,
# in this demonstration we will illustrate the principles of slice selection.
# First, let's import a 3D phantom, in this case a brain slab
# (thickness of ``2\,\mathrm{cm}``), by calling the function [`brain_phantom3D`](@ref).

obj = brain_phantom3D()
obj.Δw .= 0 # Removes the off-resonance
p1 = plot_phantom_map(obj, :T2 ; height=400)
#md savefig(p1, "../assets/3-phantom.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/3-phantom.html" style="width:50%; height:420px;"></object></center>
#md # ```

# Now, we are going to import a sequence which acquires
# 3 slices in the longitudinal axis. Note that the sequence
# contains three EPIs to acquire 3 slices of the phantom.

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/epi_multislice.seq")
seq = read_seq(seq_file)
p2 = plot_seq(seq; range=[0,10], height=400)
#md savefig(p2, "../assets/3-seq.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-seq.html" style="width:100%; height:420px;"></object>
#md # ```

# We can take a look to the slice profiles by using the function [`simulate_slice_profile`](@ref):

z = range(-2., 2., 200) * 1e-2; # -2 to 2 cm
rf1, rf2, rf3 = findall(is_RF_on.(seq))
M1 = simulate_slice_profile(seq[rf1]; z)
M2 = simulate_slice_profile(seq[rf2]; z)
M3 = simulate_slice_profile(seq[rf3]; z)

using PlotlyJS # hide
pa = scatter(x=z*1e2, y=abs.(M1.xy), name="Slice 1") # hide
pb = scatter(x=z*1e2, y=abs.(M2.xy), name="Slice 2") # hide
pc = scatter(x=z*1e2, y=abs.(M3.xy), name="Slice 3") # hide
pd = plot([pa,pb,pc], Layout(xaxis=attr(title="z [cm]"), height=300,margin=attr(t=40,l=0,r=0), title="Slice profiles for the slice-selective sequence")) # hide
#md savefig(pd, "../assets/3-profile.html") # hide
#jl display(pd)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-profile.html" style="width:100%; height:320px;"></object>
#md # ```

# Now let's simulate the acquisition.
# Notice the three echoes, one for every slice excitation.

raw = simulate(obj, seq, sys; sim_params=Dict{String,Any}("Nblocks"=>20))
p3 = plot_signal(raw; slider=false, height=300)
#md savefig(p3, "../assets/3-signal.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-signal.html" style="width:100%; height:320px;"></object>
#md # ```

# Finally, we reconstruct the acquiered images.

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
#md savefig(p4, "../assets/3-recon1.html") # hide
#md savefig(p5, "../assets/3-recon2.html") # hide
#md savefig(p6, "../assets/3-recon3.html") # hide
#jl display(p4)
#jl display(p5)
#jl display(p6)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/3-recon1.html" style="width:50%; height:380px;"></object><object type="text/html" data="../../assets/3-recon2.html" style="width:50%; height:380px;"></object>
#md # ```
#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/3-recon3.html" style="width:50%; height:380px;"></object></center>
#md # ```

# # Simple Motion Definition and Simulation

using KomaMRI # hide
sys = Scanner() # hide

# This tutorial illustrates how we can add simple motion to phantoms. We will also see how phantoms can be stored and loaded from ``.phantom`` files.

# First, we load a static 3D brain phantom:
obj = brain_phantom3D()
obj.Δw .= 0 # Removes the off-resonance

# Now, we will add Rotation Motion to recreate the patient's movement inside the scanner.

#md # !!! note
#md #     Note how rotations are defined with respect to the 3 axes:
#md #     ```@raw html
#md #     <center><img src="../../assets/head_rotation_axis.png" width="200"></center>
#md #     ```

obj.motion = SimpleMotion([Rotation(t_start=0.0, t_end=0.5, pitch=15.0, roll=0.0, yaw=45.0)])
p1 = plot_phantom_map(obj2, :T2 ; height=600, intermediate_time_samples=4)
#md savefig(p1, "../assets/5-phantom.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-phantom.html" style="width:50%; height:420px;"></object></center>
#md # ```

# Then, we will load an EPI sequence

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = read_seq(seq_file)
p2 = plot_seq(seq; range=[0 40], slider=true, height=300)
#md savefig(p3, "../assets/5-seq.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <object type="text/html" data="../../assets/5-seq.html" style="width:100%; height:320px;"></object>
#md # ```

# Now, we will run two simulations: the first with the sequence starting at ``t=0.0``, 
# and the second adding a 0.5s initial delay to the sequence:
## Simulate 
raw1 = simulate(obj, seq, sys)
raw2 = simulate(obj, Delay(0.5) + seq, sys)

# Let's note the effect of motion in both reconstructions:s
## Get the acquisition data
acq1 = AcquisitionData(raw1)
acq2 = AcquisitionData(raw2)
acq1.traj[1].circular = false #This is to remove the circular mask
acq2.traj[1].circular = false #This is to remove the circular mask

## Setting up the reconstruction parameters
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))

image1 = reconstruction(acq1, reconParams)
image2 = reconstruction(acq1, reconParams)

## Plotting the recon
p3 = plot_image(abs.(image1[:, :, 1]); height=400)
p4 = plot_image(abs.(image2[:, :, 1]); height=400)
#md savefig(p3, "../assets/5-recon1.html") # hide
#md savefig(p4, "../assets/5-recon2.html") # hide
#jl display(p3)
#jl display(p4)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-recon1.html" style="width:50%; height:420px;"></object><object type="text/html" data="../../assets/5-recon2.html" style="width:50%; height:420px;"></object></center>
#md # ```
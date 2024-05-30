# # Patient's Motion During Acquisition

using KomaMRI # hide
sys = Scanner() # hide

# It can also be interesting to see the effect of the patient's motion during an MRI scan.
# For this, Koma provides the ability to add `motion <: MotionModel` to the phantom.
# In this tutorial, we will show how to add a [`SimpleMotion`](@ref) model to a 2D brain phantom.

# First, let's load the 2D brain phantom used in the previous tutorials:
obj = brain_phantom2D()
obj.Δw .= 0 # hide

# ### Head Rotation
#
# The `SimpleMotion` model includes a list of `SimpleMotionType`'s, to enabling mix-and-matching simple motions.
# In this example, we will add a [`Rotation`](@ref) of 45 degrees around the z-axis with duration of 200 ms:

obj.motion = SimpleMotion([
    Rotation(t_start=0.0, t_end=200e-3, yaw=45.0, pitch=0.0, roll=0.0)
])
p1 = plot_phantom_map(obj, :T2 ; height=450, intermediate_time_samples=4) # hide

#md savefig(p1, "../assets/5-phantom1.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-phantom1.html" style="width:85%; height:470px;"></object></center>
#md # ```

## Read Sequence # hide
seq_file1 = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq") # hide
seq1 = read_seq(seq_file1) # hide

## Simulate # hide
raw1 = simulate(obj, seq1, sys) # hide

## Recon # hide
acq1 = AcquisitionData(raw1) # hide
acq1.traj[1].circular = false # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

# If we simulate an EPI sequence with acquisition duration (183.989 ms) comparable with the motion's duration (200 ms),
# we will observe motion-induced artifacts in the reconstructed image.
## Plotting the recon # hide
p2 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
#md savefig(p2, "../assets/5-recon1.html") # hide
#jl display(p2)

#md # ```@raw html
#md # <center>
#md # <object type="text/html" data="../../assets/5-recon1.html" style="width:65%; height:420px;"></object>
#md # </center>
#md # ```

# The severity of the artifacts can vary depending on the acquisition duration and $k$-space trajectory.

# ### Head Translation
#
# Now, let's redefine the phantom's motion with a [`Translation`](@ref) of 2 cm in x, with duration of 200 ms (v = 0.1 m/s):
obj.motion = SimpleMotion([
    Translation(t_start=0.0, t_end=200e-3, dx=2e-2, dy=0.0, dz=0.0)
])
p3 = plot_phantom_map(obj, :T2 ; height=450, intermediate_time_samples=4) # hide
#md savefig(p3, "../assets/5-phantom2.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-phantom2.html" style="width:85%; height:470px;"></object></center>
#md # ```

## Simulate # hide
raw1 = simulate(obj, seq1, sys) # hide

## Recon # hide
acq1 = AcquisitionData(raw1) # hide
acq1.traj[1].circular = false # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

# ### Motion-Corrected Reconstruction
# 
# Once simulation is done, it is possible to perform a corrected reconstrution 
# in order to revert the motion effect in the final image. 
# This can be achieved by multiplying each sample of the acquired signal 
# by a phase which is proportional to the displacement in each direction (Δx, Δy, Δz)
# at the time instant when the sample was acquired [[Godenschweger, 2016]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930872/):

# ```math
# S(k_x, k_y, k_z)_{\text{cor}} = S(k_x, k_y, k_z)_{\text{orig}} \cdot e^{i \Delta \phi_{\text{cor}}} = S(k_x, k_y, k_z)_{\text{orig}} \cdot e^{i 2 \pi (k_x \Delta x + k_y \Delta y + k_z \Delta z)}
# ```

# We need to obtain the displacements in every ADC sampling time of the sequence.
# Since translation is a rigid motion, 
# we can obtain the displacements only for one spin, 
# as the displacements of the rest will be the same.
sample_times = get_adc_sampling_times(seq1)
displacements = hcat(get_spin_coords(obj.motion, [0.0], [0.0], [0.0], sample_times)...)

p4 = KomaMRIPlots.plot( # hide
    sample_times, # hide
    displacements .* 1e2, # hide
    KomaMRIPlots.Layout( # hide
        title = "Head displacement in x, y and z", # hide
        xaxis_title = "time (s)", # hide
        yaxis_title = "Displacement (cm)" # hide
    )) # hide
KomaMRIPlots.restyle!(p4,1:3, name=["Δx", "Δy", "Δz"]) # hide

#md savefig(p4, "../assets/5-displacements.html") # hide
#jl display(p4)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-displacements.html" style="width:90%; height:470px;"></object></center>
#md # ```

# We can now get the necessary phase shift for each sample:
_, kspace = get_kspace(seq1)
ΔΦ = 2π*sum(kspace .* displacements, dims=2)

# And we apply the phase correction:
acq1.kdata[1] .*= exp.(im*ΔΦ)

image2 = reconstruction(acq1, reconParams) # hide

p5 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
p6 = plot_image(abs.(image2[:, :, 1]); height=400) # hide

#md savefig(p5, "../assets/5-recon2.html") # hide
#md savefig(p6, "../assets/5-recon3.html") # hide
#jl display(p5)
#jl display(p6)

# On the left, you can see the original reconstructed image 
# and the artifact produced by the translation in x.
# On the right, the result of the motion-corrected reconstruction, 
# where we have achieved an image similar to the one 
# we would have obtained from simulating over a static phantom.

#md # ```@raw html
#md # <object type="text/html" data="../../assets/5-recon2.html" style="width:50%; height:420px;"></object><object type="text/html" data="../../assets/5-recon3.html" style="width:50%; height:420px;"></object>
#md # ```

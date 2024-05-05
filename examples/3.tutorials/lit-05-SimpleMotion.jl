# # Patient's Motion During Acquisition

using KomaMRI # hide
sys = Scanner() # hide

# It can also be interesting to see the effect of the patient's motion during an MRI scan.
# For this, Koma provides the ability to add [`MotionModel`](@ref)'s to the phantom.
# In this tutorial, we will show how to add a [`SimpleMotion`](@ref) model to a 2D brain phantom.

# First, let's load the 2D brain phantom used in the previous tutorials:
obj = brain_phantom2D()
obj.Δw .= 0 # hide

# The `SimpleMotion` model includes a list of `SimpleMotionType`'s, to enabling mix-and-matching simple motions.
# In this example, we will add a [`Rotation`](@ref) of 20 degrees around the z-axis with duration of 200 ms.

obj.motion = SimpleMotion([
    Rotation(t_start=0.0, t_end=200e-3, yaw=20.0, pitch=0.0, roll=0.0)
    ])
p1 = plot_phantom_map(obj, :T2 ; height=400, intermediate_time_samples=4)
#md savefig(p1, "../assets/5-phantom.html") # hide
#jl display(p1)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-phantom.html" style="width:50%; height:420px;"></object></center>
#md # ```

## Read Sequence # hide
seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq") # hide
seq = read_seq(seq_file) # hide

## Simulate # hide
raw1 = simulate(obj, seq, sys) # hide

## Recon # hide
acq1 = AcquisitionData(raw1) # hide
acq1.traj[1].circular = false # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

# If we simulate an EPI sequence with acquisition duration (183.989 ms) comparable with the motion's duration (200 ms),
# we will observe motion-induced artifacts in the reconstructed image.
## Plotting the recon
p3 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
#md savefig(p3, "../assets/5-recon1.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <center>
#md # <object type="text/html" data="../../assets/5-recon1.html" style="width:50%; height:400px;"></object>
#md # </center>
#md # ```

# The severity of the artifacts can vary depending on the used acquisition duration and `k`-space trajectory.
# Below, we show the effect of the same motion in an spiral acquisition (dur. 39 ms, which is 5 times faster than the motion.)

## Read Sequence # hide
seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/5.koma_paper/comparison_accuracy/sequences/Spiral/spiral_100x100_FOV230_SPZ_INTER1.seq") # hide
seq = read_seq(seq_file) # hide

## Simulate # hide
raw1 = simulate(obj, seq, sys) # hide

## Recon # hide
acq1 = AcquisitionData(raw1) # hide
Nx, Ny = raw1.params["reconSize"][1:2] # hide
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny)) # hide
image1 = reconstruction(acq1, reconParams) # hide

## Plotting the recon # hide
p4 = plot_image(abs.(image1[:, :, 1]); height=400) # hide
#md savefig(p4, "../assets/5-recon2.html") # hide
#jl display(p4)

#md # ```@raw html
#md # <center>
#md # <object type="text/html" data="../../assets/5-recon2.html" style="width:50%; height:400px;"></object>
#md # </center>
#md # ```

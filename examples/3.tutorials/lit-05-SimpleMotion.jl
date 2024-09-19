# # Patient's Motion During Acquisition

using KomaMRI # hide
sys = Scanner() # hide

# It can also be interesting to see the effect of the patient's motion during an MRI scan.
# For this, Koma provides the ability to add `motion <: AbstractMotion` to the phantom.
# In this tutorial, we will show how to add a [`Translate`](@ref) motion to a 2D brain phantom.

# First, let's load the 2D brain phantom used in the previous tutorials:
obj = brain_phantom2D()
obj.Δw .= 0 # hide

# ### Head Translation
#
# In this example, we will add a [`Translate`](@ref) of 2 cm in x, with duration of 200 ms (v = 0.1 m/s):

obj.motion = MotionList(
    Translate(2e-2, 0.0, 0.0, TimeRange(t_start=0.0, t_end=200e-3))
)
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

# ### Motion-Corrected Reconstruction
# 
# To correct for the motion-induced artifacts we can perform a motion-corrected reconstruction. 
# This can be achieved by multiplying each sample of the acquired signal  ``S(t)``
# by a phase shift ``\Delta\phi_{\mathrm{corr}}`` proportional to the displacement ``\boldsymbol{u}(t)``
# [[Godenschweger, 2016]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930872/):

# ```math
# S_{\mathrm{MC}}\left(t\right)=S\left(t\right)\cdot\mathrm{e}^{\mathrm{i}\Delta\phi_{\mathrm{corr}}}=S\left(t\right)\cdot\mathrm{e}^{\mathrm{i}2\pi\boldsymbol{k}\left(t\right)\cdot\boldsymbol{u}\left(t\right)}
# ```

# In practice, we would need to estimate or measure the motion before performing a motion-corrected reconstruction, but for this example, we will directly use the displacement functions ``\boldsymbol{u}(\boldsymbol{x}, t)`` defined by `obj.motion::MotionList`. 
# Since translations are rigid motions (``\boldsymbol{u}(\boldsymbol{x}, t)=\boldsymbol{u}(t)`` no position dependence), we can obtain the required displacements by calculating ``\boldsymbol{u}(\boldsymbol{x}=\boldsymbol{0},\ t=t_{\mathrm{adc}})``.
sample_times = get_adc_sampling_times(seq1)
displacements = hcat(get_spin_coords(obj.motion, [0.0], [0.0], [0.0], sample_times)...)

p3 = KomaMRIPlots.plot( # hide
    sample_times, # hide
    displacements .* 1e2, # hide
    KomaMRIPlots.Layout( # hide
        title = "Head displacement in x, y and z", # hide
        xaxis_title = "time (s)", # hide
        yaxis_title = "Displacement (cm)" # hide
    )) # hide
KomaMRIPlots.restyle!(p3,1:3, name=["ux(t)", "uy(t)", "uz(t)"]) # hide

#md savefig(p3, "../assets/5-displacements.html") # hide
#jl display(p3)

#md # ```@raw html
#md # <center><object type="text/html" data="../../assets/5-displacements.html" style="width:80%; height:300px;"></object></center>
#md # ```

# We can now get the necessary phase shift for each sample:
_, kspace = get_kspace(seq1)
ΔΦ = 2π*sum(kspace .* displacements, dims=2)

# And apply it to the acquired signal to correct its phase:
acq1.kdata[1] .*= exp.(im*ΔΦ)

image2 = reconstruction(acq1, reconParams) # hide

p4 = plot_image(abs.(image2[:, :, 1]); height=400) # hide

#md savefig(p4, "../assets/5-recon2.html") # hide

#jl display(p2)
#jl display(p4)

# Finally, we compare the original image ▶️ and the motion-corrected reconstruction ⏸️:

#md # ```@raw html
#md # <object type="text/html" data="../../assets/5-recon1.html" style="width:50%; height:420px;"></object><object type="text/html" data="../../assets/5-recon2.html" style="width:50%; height:420px;"></object>
#md # ```

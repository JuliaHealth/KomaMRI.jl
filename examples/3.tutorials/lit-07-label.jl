# # Using Labels to reconstruct multi-slice / multi-constrast sequences
# 
# This is generally done in pulseq / pypulseq. In that case the sequence structure automatically contains the labels increment and set extension.
#
# Here we will show how to modify a sequence in KomaMRI to add the extension

using KomaMRI, PlotlyJS # hide
sys = Scanner() # hide
obj = brain_phantom3D()

# ## How to add label to a sequence
# let's add the label to increment the slice for each seq_EPI block. We will add that to the first block of the sequence

seq_EPI = PulseDesigner.EPI_example()
seq = deepcopy(seq_EPI)
# We can create 2 different types of label "Increment" and "Set" with `LabelInc` and `LabelSet`
lSlcInc = LabelInc(1,"SLC")

# Let's change the extension field of the first block of seq_EPI in order to add an increment to the SLC label
seq_EPI.EXT[1] = [lSlcInc]

## now let's merge 3 seq_EPI we now have 1 EPI sequence without EXTENSION then we have an increment of the SLC label and another one at the begining of the last seq_EPI
seq = seq + seq_EPI + seq_EPI
plot_seq(seq)
# we can extract the label value with the following fonction
l = get_label(seq)
# and store in a vector only the value of the SLC label
SLC_vec = [l[i].SLC for i in eachindex(l)] 

# The value is equal to 0 until we reach LabelInc(1,"SLC) 
plot(SLC_vec, Layout(
    xaxis_title="n° blocks",
    yaxis_title="SLC label"
))


# ## Simulate the acquisition and reconstruct the data
# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
raw = simulate(obj, seq, sys; sim_params)

# the simulated rawdata stored the correct label for each profiles 
for p in [1,101,102,203]
  slc = raw.profiles[p].head.idx.slice |> Int 
  println("Profile n°$p : SLC label = $slc")
end

#  MRIReco splits the data in 3 when the data are converted to the AcquisitionData structure. The dimension are (contrast / slice / repetitions )
acqData = AcquisitionData(raw)
acqData.kdata |> size

# For multi-slice acquisition MRIReco uses the same trajectory. We need to crop the trajectory.
size_readout = size(acqData.traj[1].nodes,2) / 3 |> Int
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,1:size_readout]
# We need to normalize the trajectory to -0.5 to 0.5
C = maximum(2*abs.(acqData.traj[1].nodes[1:2]))
acqData.traj[1].nodes ./= C
acqData.traj[1].circular = false #Removing circular window

Nx, Ny = raw.params["reconSize"][1:2]

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
rec = reconstruction(acqData, recParams)
image=abs.(reshape(rec,Nx,Ny,:))

# Let's take a look at the 3 images

plot_image(image[:,:,1])
plot_image(image[:,:,2])
plot_image(image[:,:,3])

# The signal ponderation is changing because we are acquired the same slice position with a short TR sequence. Thus, the magnetization is not at the equilibrium.

# ## 
## Rather than using the k-space trajectory calculated by KomaMRI and perform a NUFFT for the reconstruction, we can use the label LIN (phase encoding) and PAR (partition encoding).
#
# First we will create an increment label for LIN
seq_LIN =  PulseDesigner.EPI_example()
lInc = LabelInc(1,"LIN")

# We can put the increment at anystage between 2 ADC blocks. Here we will put it onto each ADC blocks.
idx_ADC = is_ADC_on.(seq_LIN)
for i in eachindex(idx_ADC)
  idx_ADC[i] == 1 ? seq_LIN.EXT[i] = [lInc] : nothing
end

# because we want the label of each ADC starting from 0. We set the value to -1 on the first block.
seq_LIN.EXT[1] = [LabelSet(-1,"LIN")]

# let's check the LIN label for each adc
l = get_label(seq_LIN)
l[idx_ADC]

# we can now simulate the results
raw = simulate(obj, seq_LIN, sys; sim_params)

# In order to not use the NUFFT reconstruction of MRIReco.jl we need to change the trajectory name to "cartesian"
raw.params["trajectory"] = "cartesian"
raw.params["encodedSize"] = [seq.DEF["Nx"],seq.DEF["Ny"]]

# We also need to estimate the profile center which will be at the center of the readout. If it is not the case it should be specified in `f.profiles[i].head.center_sample = center_sample` and `estimateProfileCenter = false`
acqData = AcquisitionData(raw,estimateProfileCenter=true)

# Let's see the results
recParams = Dict{Symbol,Any}()
rec = reconstruction(acqData, recParams)

plot_image(abs.(rec[:,:,1]))

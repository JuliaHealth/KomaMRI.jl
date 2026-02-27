using KomaMRI, Suppressor #hide
sys = Scanner(); #hide

obj = brain_phantom3D()
obj.Δw .= 0 # Removes the off-resonance
p1 = @suppress plot_phantom_map(obj, :T2 ; height=400)
display(p1);

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/epi_multislice.seq")
seq = @suppress read_seq(seq_file)
p2 = plot_seq(seq; range=[0,10], height=400)
display(p2);

z = range(-2., 2., 200) * 1e-2; # -2 to 2 cm
rf1, rf2, rf3 = findall(is_RF_on.(seq))
M1 = @suppress simulate_slice_profile(seq[rf1]; z)
M2 = @suppress simulate_slice_profile(seq[rf2]; z)
M3 = @suppress simulate_slice_profile(seq[rf3]; z);

using PlotlyJS #hide
pa = scatter(x=z*1e2, y=abs.(M1.xy), name="Slice 1") #hide
pb = scatter(x=z*1e2, y=abs.(M2.xy), name="Slice 2") #hide
pc = scatter(x=z*1e2, y=abs.(M3.xy), name="Slice 3") #hide
pd = plot([pa,pb,pc], Layout(xaxis=attr(title="z [cm]"), height=300,margin=attr(t=40,l=0,r=0), title="Slice profiles for the slice-selective sequence")) #hide
display(pd)

raw = @suppress simulate(obj, seq, sys; sim_params=Dict{String,Any}("Nblocks"=>20))
p3 = plot_signal(raw; slider=false, height=300)
display(p3)

# Get the acquisition data
acq = AcquisitionData(raw)

# Setting up the reconstruction parameters and perform reconstruction
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
image = reconstruction(acq, reconParams)

# Plotting the slices
p4 = plot_image(abs.(image[:, :, 1]); height=360, title="Slice 1")
p5 = plot_image(abs.(image[:, :, 2]); height=360, title="Slice 2")
p6 = plot_image(abs.(image[:, :, 3]); height=360, title="Slice 3")
display([p4 p5 p6])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

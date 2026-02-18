using KomaMRI, PlotlyJS, Suppressor # hide
sys = Scanner() # hide
obj = brain_phantom3D()
obj.Δw .= 0; # Removes the off-resonance

seq_EPI = PulseDesigner.EPI_example()
seq = deepcopy(seq_EPI);

lSlcInc = LabelInc(1,"SLC");

seq_EPI.EXT[1] = [lSlcInc];

seq = seq + seq_EPI + seq_EPI
plot_seq(seq);

l = get_label(seq);

SLC_vec = [l[i].SLC for i in eachindex(l)];

p1 = plot(SLC_vec, Layout(
    xaxis_title="n° blocks",
    yaxis_title="SLC label"
))
display(p1);

sim_params = KomaMRICore.default_sim_params()
raw = @suppress simulate(obj, seq, sys; sim_params)

for p in [1,101,102,203]
  slc = raw.profiles[p].head.idx.slice |> Int
  println("Profile n°$p : SLC label = $slc")
end

acqData = AcquisitionData(raw)
acqData.kdata |> size

size_readout = size(acqData.traj[1].nodes,2) / 3 |> Int
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,1:size_readout];

C = maximum(2*abs.(acqData.traj[1].nodes[1:2]))
acqData.traj[1].nodes ./= C
acqData.traj[1].circular = false # Removing circular window

Nx, Ny = raw.params["reconSize"][1:2]

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
rec = reconstruction(acqData, recParams)
image = abs.(reshape(rec,Nx,Ny,:));

p2 = plot_image(image[:,:,1], height=400);
p3 = plot_image(image[:,:,2], height=400);
display(p2);
display(p3);

seq_LIN = PulseDesigner.EPI_example()
lInc = LabelInc(1,"LIN");

idx_ADC = is_ADC_on.(seq_LIN)
for i in eachindex(idx_ADC)
  idx_ADC[i] == 1 ? seq_LIN.EXT[i] = [lInc] : nothing;
end

seq_LIN.EXT[1] = [LabelSet(-1,"LIN")];

l = get_label(seq_LIN)
l[idx_ADC][1:10]

raw = @suppress simulate(obj, seq_LIN, sys; sim_params)

raw.params["trajectory"] = "cartesian"
raw.params["encodedSize"] = [seq.DEF["Nx"],seq.DEF["Ny"]];

acqData = AcquisitionData(raw,estimateProfileCenter=true);

recParams = Dict{Symbol,Any}()
rec = reconstruction(acqData, recParams);

p4=plot_image(abs.(rec[:,:,1]), height=400);
display(p4);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

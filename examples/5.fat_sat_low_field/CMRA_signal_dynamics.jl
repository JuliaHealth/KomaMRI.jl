using  KomaMRI, PlotlyJS

# Phantom
include("CMRA_phantom.jl")
# Scanner
sys = Scanner()
# Sequence
include("CMRA_seq.jl")
# Simulation
seq = CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines,
                        im_segments, iNAV_flip_angle, im_flip_angle,
                        T2prep_duration, FatSat_flip_angle;
                        sample_recovery=ones(Bool,number_dummy_heart_beats+1))
obj = cardiac_phantom(0)
sim_method = BlochDict(save_Mz=true)
sim_params = Dict{String,Any}("return_type"=>"mat", "sim_method"=>sim_method, "gpu"=>false,
                             "Nthreads"=>1, "Δt_rf"=>2e-4)
magnetization = simulate(obj, seq, sys; sim_params)
# Prep plots
labels = ["Myocardium", "Blood", "Fat"]
colors = ["blue", "red", "green"]
spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)']
t = KomaMRICore.get_adc_sampling_times(seq)
Mxy(i) = abs.(sum(magnetization[:,spin_group[i],1,1][:,1,:],dims=2)[:]/length(spin_group[i]))
Mz(i) = real.(sum(magnetization[:,spin_group[i],2,1][:,1,:],dims=2)[:]/length(spin_group[i]))
# Plot
p0 = make_subplots(rows=3, cols=1, subplot_titles=["M_xy" "Mz" "Sequence"], shared_xaxes=true, vertical_spacing=0.1)
for i=eachindex(spin_group)
    p1 = scatter(x=t, y=Mxy(i), name=labels[i], legendgroup=labels[i], marker_color=colors[i])
    p2 = scatter(x=t, y=Mz(i), name=labels[i], legendgroup=labels[i], showlegend=false, marker_color=colors[i])#,line=attr(dash="dash"))
    add_trace!(p0, p1, row=1, col=1)
    add_trace!(p0, p2, row=2, col=1)
end
seqd = KomaMRICore.discretize(seq; sampling_params=sim_params)
p3 = scatter(x=seqd.t, y=abs.(seqd.B1), name="B1",marker_color="purple",yaxis_range=[0,5])
add_trace!(p0, p3, row=3, col=1)
add_layout_image!(p0, attr(
                    source="https://raw.githubusercontent.com/cncastillo/KomaMRI.jl/master/assets/logo.svg",
                    xref="x domain",
                    yref="y domain",
                    x=0.99,
                    y=1.2,
                    opacity=0.7,
                    xanchor="right",
                    yanchor="top",
                    sizex=0.15,
                    sizey=0.15,))
relayout!(p0, yaxis_range=[0, 0.4],
        xaxis_range=[RR*number_dummy_heart_beats, RR*number_dummy_heart_beats+.250],
        title_text="TR=$(round(TR*1e3;digits=3)) ms, α=$(im_flip_angle), iNAV_lines=$(iNAV_lines), FatSat α=$(FatSat_flip_angle)")
p0

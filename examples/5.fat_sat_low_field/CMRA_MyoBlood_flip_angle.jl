using  KomaMRI, PlotlyJS

# Phantom
include("CMRA_phantom.jl")
# Scanner
sys = Scanner()
# Sequence
include("CMRA_seq.jl")
# Simulation
RRs = 0.7 : 0.05 : 1.1 #s
param = (20:5:180) #* 1e-3
magnetization = zeros(ComplexF64, im_segments, Niso*3, length(param), length(RRs))
for (m, RR) = enumerate(RRs)
    for (n, im_flip_angle) = enumerate(param)
        seq = CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, 
                                im_segments, iNAV_flip_angle, im_flip_angle,
                                T2prep_duration, FatSat_flip_angle, RR; 
                                sample_recovery=zeros(Bool,number_dummy_heart_beats+1))
        obj = cardiac_phantom(0)
        sim_method = BlochDict(save_Mz=false)
        simParams = Dict{String,Any}("return_type"=>"mat", "sim_method"=>sim_method, "gpu"=>false, 
                                    "Nthreads"=>1, "Δt_rf"=>2e-4)
        magnetization[:, :, n, m] = simulate(obj, seq, sys; simParams)[end-im_segments+1:end, :]
    end
end
## Plotting results
labels = ["Myocardium", "Blood", "Fat"]
colors = ["blue", "red", "green", "purple"]
spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)'] 
mean(x, dim) = sum(x, dims=dim) / size(x, dim)
std(x, dim; mu=mean(x,dim)) = sqrt.(sum(abs.(x .- mu).^2, dims=dim) / (size(x,dim) - 1))
signal_myoc = reshape(mean(abs.(mean(magnetization[:, spin_group[1], :, :], 3)), 1), length(param), length(RRs))
signal_bloo = reshape(mean(abs.(mean(magnetization[:, spin_group[2], :, :], 3)), 1), length(param), length(RRs))
signal_fat  = reshape(mean(abs.(mean(magnetization[:, spin_group[3], :, :], 3)), 1), length(param), length(RRs))
diff_bloo_myoc = abs.(signal_bloo .- signal_myoc)
mean_myoc = mean(signal_myoc,2)
mean_bloo = mean(signal_bloo,2)
mean_fat  = mean(signal_fat, 2)
mean_diff = mean(diff_bloo_myoc,2)
std_myoc = std(signal_myoc,2) 
std_bloo = std(signal_bloo,2)
std_fat  = std(signal_fat, 2) 
std_diff = std(diff_bloo_myoc,2)
# Plotting results
~, idx1 = findmax(mean_bloo); param[idx1]
~, idx2 = findmax(mean_myoc); param[idx2]
~, idx3 = findmax(mean_diff); param[idx3]
s1 = scatter(x=param, 
             y=mean_myoc[:],
             name=labels[1],legendgroup=labels[1], line=attr(color=colors[1]))
s2 = scatter(x=param, 
             y=mean_bloo[:],
             name=labels[2],legendgroup=labels[2], line=attr(color=colors[2]))
s3 = scatter(x=param, 
             y=mean_fat[:],
             name=labels[3],legendgroup=labels[3], line=attr(color=colors[3]))
s4 = scatter(x=param, 
             y=mean_diff[:],
             name="|Blood-Myoc|",legendgroup="|Blood-Myoc|", line=attr(color=colors[4]))
s5 = scatter(x=[param; reverse(param)], 
             y=[(mean_myoc .- std_myoc)[:]; reverse((mean_myoc .+ std_myoc)[:])],
             name=labels[1],legendgroup=labels[1],showlegend=false,fill="toself",
             fillcolor="rgba(0,0,255,0.2)", line=attr(color="rgba(0,0,0,0)"))
s6 = scatter(x=[param; reverse(param)], 
             y=[(mean_bloo .- std_bloo)[:]; reverse((mean_bloo .+ std_bloo)[:])],
             name=labels[2],legendgroup=labels[2],showlegend=false,fill="toself",
             fillcolor="rgba(255,0,0,0.2)", line=attr(color="rgba(0,0,0,0)"))
s7 = scatter(x=[param; reverse(param)], 
             y=[(mean_fat .- std_fat)[:]; reverse((mean_fat .+ std_fat)[:])],
             name=labels[3],legendgroup=labels[3],showlegend=false,fill="toself",
             fillcolor="rgba(0,255,0,0.2)", line=attr(color="rgba(0,0,0,0)"))
s8 = scatter(x=[param; reverse(param)], 
             y=[(mean_diff .- std_diff)[:]; reverse((mean_diff .+ std_diff)[:])],
             name="|Blood-Myoc|",legendgroup="|Blood-Myoc|",showlegend=false,fill="toself",
             fillcolor="rgba(255,0,255,0.2)", line=attr(color="rgba(0,0,0,0)"))
fig = plot([s1,s2,s3,s4,s5,s6,s7,s8])
relayout!(fig, yaxis=attr(title="Signal [a.u.]",tickmode="array", 
               tickvals=round.([0,mean_myoc[idx2],mean_bloo[idx1],mean_diff[idx3]], digits=3),
               ticktext=round.([0,mean_myoc[idx2],mean_bloo[idx1],mean_diff[idx3]], digits=2)
               ),
               font=attr(family="CMU Serif", size=16, 
               scaleanchor="x", scaleratio=1),
               xaxis=attr(title="Flip angle [deg]",
               tickmode="array", 
               tickvals=[param[idx2],param[idx1],param[idx3],param[end]],
               constrain="domain"),
               yaxis_range=[0,.3],
               width=500, height=400,
               hovermode = "x unified")
savefig(fig, "/home/ccp/Desktop/CMRA_FA.svg")
fig
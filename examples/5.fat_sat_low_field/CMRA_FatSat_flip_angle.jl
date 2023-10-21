using  KomaMRI, PlotlyJS

# Phantom
include("CMRA_phantom.jl")
# Scanner
sys = Scanner()
# Sequence
include("CMRA_seq.jl")
# Simulation
inhom_ppm = 1e-6 # From paper 0.5 ppm
inhom_freq_shift = γ * B0 * inhom_ppm
offs = -inhom_freq_shift:inhom_freq_shift  #s
param = 20:10:250 #flip angle [deg]
magnetization = zeros(ComplexF64, im_segments, Niso*4, length(param), length(offs))
for (m, off) = enumerate(offs)
    for (n, FatSat_flip_angle) = enumerate(param)
        seq = CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, 
                                im_segments, iNAV_flip_angle, im_flip_angle,
                                T2prep_duration, FatSat_flip_angle, RR; 
                                sample_recovery=zeros(Bool,number_dummy_heart_beats+1))
        obj = cardiac_phantom(off)
        simParams = Dict{String,Any}(
            "return_type"=>"mat", 
            "sim_method"=>BlochDict(save_Mz=false), 
            "gpu"=>false, 
            "Nthreads"=>1, 
            "Δt_rf"=>1.7e-4
            )
        magnetization[:, :, n, m] = simulate(obj, seq, sys; simParams)[end-im_segments+1:end, :]
    end
end
## Calculating results
labels = ["Myocardium", "Blood", "Fat (T₁=183 ms)", "Fat (T₁=130 ms)"]
colors = ["blue", "red", "green"]
spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)', (3Niso+1:4Niso)'] 
mean(x, dim) = sum(x, dims=dim) / size(x, dim)
std(x, dim; mu=mean(x, dim)) = sqrt.(sum(abs.(x .- mu).^2, dims=dim) / (size(x, dim) - 1))
signal_myoc = reshape(mean(abs.(mean(magnetization[:, spin_group[1], :, :], 3)), 1), length(param), length(offs))
signal_bloo = reshape(mean(abs.(mean(magnetization[:, spin_group[2], :, :], 3)), 1), length(param), length(offs))
signal_fat  = reshape(mean(abs.(mean(magnetization[:, spin_group[3], :, :], 3)), 1), length(param), length(offs))
signal_fat2 = reshape(mean(abs.(mean(magnetization[:, spin_group[4], :, :], 3)), 1), length(param), length(offs))
mean_myoc = mean(signal_myoc,2)
mean_bloo = mean(signal_bloo,2)
mean_fat  = mean(signal_fat, 2)
mean_fat2  = mean(signal_fat2, 2)
std_myoc =  std(signal_myoc,2) 
std_bloo =  std(signal_bloo,2) 
std_fat  =  std(signal_fat, 2) 
std_fat2 =  std(signal_fat2, 2) 
# Plotting results
~, idx1 = findmax(mean_bloo); param[idx1]
~, idx2 = findmax(mean_myoc); param[idx2]
~, idx3 = findmax(mean_fat ); param[idx3]
s1 = scatter(x=param, 
             y=mean_myoc[:],
             name=labels[1],legendgroup=labels[1], line=attr(color=colors[1]))
s2 = scatter(x=param, 
             y=mean_bloo[:],
             name=labels[2],legendgroup=labels[2], line=attr(color=colors[2]))
s3 = scatter(x=param, 
             y=mean_fat[:],
             name=labels[3],legendgroup=labels[3], line=attr(color=colors[3]))
s4 = scatter(x=[param; reverse(param)], 
             y=[(mean_myoc .- std_myoc)[:]; reverse((mean_myoc .+ std_myoc)[:])],
             name=labels[1],legendgroup=labels[1],showlegend=false,fill="toself",
             fillcolor="rgba(0,0,255,0.2)", line=attr(color="rgba(0,0,0,0)"))
s5 = scatter(x=param, 
             y=mean_fat2[:],
             name=labels[4],legendgroup=labels[4], line=attr(color=colors[3], dash="dash"))
s6 = scatter(x=[param; reverse(param)], 
             y=[(mean_bloo .- std_bloo)[:]; reverse((mean_bloo .+ std_bloo)[:])],
             name=labels[2],legendgroup=labels[2],showlegend=false,fill="toself",
             fillcolor="rgba(255,0,0,0.2)", line=attr(color="rgba(0,0,0,0)"))
s7 = scatter(x=[param; reverse(param)], 
             y=[(mean_fat .- std_fat)[:]; reverse((mean_fat .+ std_fat)[:])],
             name=labels[3],legendgroup=labels[3],showlegend=false,fill="toself",
             fillcolor="rgba(0,255,0,0.2)", line=attr(color="rgba(0,0,0,0)"))
s8 = scatter(x=[param; reverse(param)], 
             y=[(mean_fat2 .- std_fat2)[:]; reverse((mean_fat2 .+ std_fat2)[:])],
             name=labels[4],legendgroup=labels[4],showlegend=false,fill="toself",
             fillcolor="rgba(0,150,100,0.2)", line=attr(color="rgba(0,0,0,0)"))

fig = plot([s1,s2,s3,s4,s5,s6,s7,s8])
relayout!(fig, 
            yaxis=attr(
                title="Signal [a.u.]",tickmode="array", 
                # tickvals=round.([0,mean_myoc[], digits=3),
                # ticktext=round.([0,mean_myoc[], digits=2)
            ),
            xaxis=attr(
                title="FatSat flip angle [deg]",
                tickmode="array", 
                tickvals=[param[1], 130, 150, 180, param[end]],
                constrain="domain"
            ),
            font=attr(
                family="CMU Serif", size=16, 
                scaleanchor="x", scaleratio=1
            ),
            yaxis_range=[0,.4],
            width=600, height=400,
            hovermode = "x unified")
fig
## Saving figure
savefig(fig, "/home/ccp/Desktop/CMRA_FA_FatSat.svg")
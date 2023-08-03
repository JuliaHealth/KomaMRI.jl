using  KomaMRI, PlotlyJS, MAT
#Phantom, values from: 
# A review of normal tissue hydrogen NMR relaxation times and relaxation mechanisms from 1–100 MHz: 
# Dependence on tissue type, NMR frequency, temperature, species, excision, and age, Paul A. Bottomley,  
# Thomas H. Foster,  Raymond E. Argersinger,  Leah M. Pfeifer
B0 = 0.55
fat_ppm = -3.4e-6
Niso = 200
Δx_voxel = 1.5e-3
fat_freq = γ*B0*fat_ppm
off = Array(range(-1, 1, Niso))
dx = Array(range(-Δx_voxel/2, Δx_voxel/2, Niso))
#BASED ON T1 and T2 maps
myocard = Phantom{Float64}(x=dx, T1=700e-3*ones(Niso), T2=60e-3*ones(Niso), Δw=2π*20*off)
blood =   Phantom{Float64}(x=dx, T1=1250e-3*ones(Niso), T2=300e-3*ones(Niso), Δw=2π*20*off) #1.5 T2 240ms
fat =     Phantom{Float64}(x=dx, T1=172.93e-3*ones(Niso), T2=135.32e-3*ones(Niso), Δw=2π*(fat_freq .+ 20*off))
obj = myocard+blood+fat
#Sequence parameters
TR = 5.29e-3
# TR = 4.89e-3
im_flip_angle = 110
iNAV_lines = 14
FatSat_flip_angle = 110
#Pre-pulses
T2prep_duration = 50e-3
Tfatsat = 32e-3
Δf_sinc = fat_freq
additional_offset = 0
# Rest of parameters
TE = TR / 2 #bSSFP condition
Tadc = 1e-6
RR = 1 #1 [s]
Trf = 1.4e-3  #1 [ms]
B1 = 1 / (360*γ*Trf)
im_segments = 20
iNAV_flip_angle = 3.2
number_dummy_heart_beats = 3
#Scanner
sys = Scanner()

function FatSat(α, Δf; sample=false)
    RF_wf = matread("./examples/5.fat_sat_low_field/GAUSS5120.mat")["B1"]
    seq = Sequence()
    seq += RF(RF_wf, Tfatsat, Δf)
    α_ref = get_flip_angles(seq)[2]
    seq *= (α/α_ref + 0im)
    if sample
        seq += ADC(1, 1e-6)
    end
    seq += Grad(60e-3, 4.7e-3, 4e-3)
    if sample
        seq += ADC(1, 1e-6)
    end
    return seq
end

function T2prep(TE; sample=false)
    seq = Sequence()
    seq += RF(90 * B1, Trf)
    seq += sample ? ADC(20, TE/2) : Delay(TE/2)
    seq += RF(180im * B1/2, 2Trf)
    seq += sample ? ADC(20, TE/2) : Delay(TE/2)
    seq += RF(-90 * B1, Trf)
    seq += Grad(30e-3, 4.7e-3, 4e-3)
    if sample
        seq += ADC(1, 1e-6)
    end
    return seq
end

function bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=false)
    k = 0
    seq = Sequence()
    for i = 0 : iNAV_lines + im_segments - 1
        if iNAV_lines != 0
            m = (im_flip_angle - iNAV_flip_angle) / iNAV_lines
            α = min( m * i + iNAV_flip_angle, im_flip_angle ) * exp(1im * π * k)
        else
            α = im_flip_angle * exp(1im * π * k)
        end
        seq += RF(α * B1, Trf)
        if i < iNAV_lines && !sample
            seq += Delay(TR - Trf)
        else
            seq += Delay(TE - Trf/2 - Tadc/2)
            seq += ADC(1, Tadc)
            seq += Delay(TR - TE - Tadc/2 - Trf/2)
        end
        k += 1
    end
    return seq
end

function CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample_recovery=ones(Bool,number_dummy_heart_beats+1))
    k = 0
    seq = Sequence()
    for hb = 0 : number_dummy_heart_beats
        t2p = T2prep(T2prep_duration; sample=sample_recovery[hb+1])
        fatsat = FatSat(FatSat_flip_angle, Δf_sinc + additional_offset; sample=sample_recovery[hb+1])
        bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=sample_recovery[hb+1])
        RRdelay = RR  - dur(bssfp) - dur(t2p) - dur(fatsat)
        seq += t2p
        seq += fatsat
        seq += bssfp
        seq += sample_recovery[hb+1] ? ADC(80, RRdelay) : Delay(RRdelay) #Sampling recovery curve
    end
    return seq
end

#Sequence
seq = CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, 
                        im_segments, iNAV_flip_angle, im_flip_angle; sample_recovery=[false, false, false, true])
#Simulation
sim_method = BlochDict(save_Mz=true)
simParams = Dict{String,Any}("return_type"=>"mat", "sim_method"=>sim_method, "gpu"=>false, "Nthreads"=>1)
magnetization = simulate(obj, seq, sys; simParams)

labels = ["Myocardium", "Blood", "Fat"]
colors = ["blue", "red", "green"]
spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)'] 
t = KomaMRICore.get_adc_sampling_times(seq) 

Mxy(i) = abs.(sum(magnetization[:,spin_group[i],1,1][:,1,:],dims=2)[:]/length(spin_group[i]))
Mz(i) = real.(sum(magnetization[:,spin_group[i],2,1][:,1,:],dims=2)[:]/length(spin_group[i]))
p0 = make_subplots(rows=3, cols=1, subplot_titles=["Mxy" "Mz" "Sequence"], shared_xaxes=true, vertical_spacing=0.1)
for i=eachindex(spin_group)
    p1 = scatter(x=t, y=Mxy(i), name=labels[i], legendgroup=labels[i], marker_color=colors[i])
    p2 = scatter(x=t, y=Mz(i), name=labels[i], legendgroup=labels[i], showlegend=false, marker_color=colors[i])#,line=attr(dash="dash"))
    add_trace!(p0, p1, row=1, col=1)
    add_trace!(p0, p2, row=2, col=1)
end
seqd = KomaMRICore.discretize(seq; simParams)
p3 = scatter(x=seqd.t, y=abs.(seqd.B1)*1e6, name="B1",marker_color="purple",yaxis_range=[0,5])
add_trace!(p0, p3, row=3, col=1)
add_layout_image!(p0, attr(
                    source="https://raw.githubusercontent.com/cncastillo/KomaMRI.jl/master/src/ui/assets/Logo.svg",
                    xref="x domain",
                    yref="y domain",
                    x=0.99,
                    y=1.2,
                    opacity=0.7,
                    xanchor="right",
                    yanchor="top",
                    sizex=0.15,
                    sizey=0.15,))
relayout!(p0, yaxis_range=[0, 0.5], xaxis_range=[3, 4],title_text="TR=$(round(TR*1e3;digits=3)) ms, α=$(im_flip_angle), iNAV_lines=$(iNAV_lines), FatSat α=$(FatSat_flip_angle)")
p0

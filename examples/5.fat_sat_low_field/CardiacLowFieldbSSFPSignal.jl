using  KomaMRI, PlotlyJS, MAT
#Phantom, values from: 
# A review of normal tissue hydrogen NMR relaxation times and relaxation mechanisms from 1–100 MHz: 
# Dependence on tissue type, NMR frequency, temperature, species, excision, and age, Paul A. Bottomley,  
# Thomas H. Foster,  Raymond E. Argersinger,  Leah M. Pfeifer
B0 = 0.55
fat_ppm = -3.4e-6
Niso = 200
Δx_voxel = 2e-3
fat_freq = γ*B0*fat_ppm
off = Array(range(-1, 1, Niso))
dx = Array(range(-Δx_voxel/2, Δx_voxel/2, Niso))
#BASED ON T1 and T2 maps
myocard = Phantom{Float64}(x=dx, T1=700e-3*ones(Niso), T2=60e-3*ones(Niso))
blood =   Phantom{Float64}(x=dx, T1=1250e-3*ones(Niso), T2=300e-3*ones(Niso)) #1.5 T2 240ms
fat =     Phantom{Float64}(x=dx, T1=96e-3*ones(Niso), T2=150e-3*ones(Niso), Δw=2π*(fat_freq .+ 5*off))
#OLD VALUES, to Ivan
# myocard = Phantom{Float64}(x=dx*0, T1=600e-3*ones(Niso), T2=60e-3*ones(Niso))
# blood =   Phantom{Float64}(x=dx*0, T1=1200e-3*ones(Niso),T2=410e-3*ones(Niso))
# fat =     Phantom{Float64}(x=dx*0, T1=217e-3*ones(Niso), T2=95e-3*ones(Niso), Δw=2π*(fat_freq .+ 5*off))
obj = myocard+blood+fat
#Sequence parameters
TR = 4.4e-3 #5.81
TE = TR / 2 #bSSFP condition
Tadc = 1e-6
RR = 1 #1 [s]
Trf = 1e-3  #1 [ms]
B1 = 1 / (360*γ*Trf)
iNAV_lines = 6
im_segments = 22
iNAV_flip_angle = 3.2
im_flip_angle = 80 #90
number_dummy_heart_beats = 3
#Pre-pulses
T2prep_duration = 50e-3
Tfatsat = 32e-3
Δf_sinc = fat_freq
additional_offset = 0
FatSat_flip_angle = 180
#Scanner
sys = Scanner()

function FatSat(α, Δf; sample=false)
    RF_wf = matread("./examples/5.fat_sat_low_field/GAUSS5120.mat")["B1"]
    seq = Sequence()
    seq += RF(RF_wf, Tfatsat, Δf)
    α_ref = get_flip_angles(seq)[2]
    seq *= (α/α_ref + 0im)  
    seq += Grad(20e-3, 4.7e-3, 1e-3)
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
    seq += Grad(20e-3, 4.7e-3, 1e-3)
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

function CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample_recovery=false)
    k = 0
    seq = Sequence()
    for hb = 0 : number_dummy_heart_beats
        t2p = T2prep(T2prep_duration; sample=sample_recovery)
        fatsat = FatSat(FatSat_flip_angle, Δf_sinc + additional_offset; sample=sample_recovery)
        bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=sample_recovery)
        RRdelay = RR  - dur(bssfp) - dur(t2p) - dur(fatsat)
        seq += t2p
        seq += fatsat             
        seq += bssfp
        seq += sample_recovery ? ADC(40, RRdelay) : Delay(RRdelay) #Sampling recovery curve
    end
    return seq
end

function BOOST_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample_recovery=false)
    k = 0
    seq = Sequence()
    for hb = 0 : number_dummy_heart_beats
        t2p = T2prep(T2prep_duration; sample=sample_recovery)
        fatsat = FatSat(FatSat_flip_angle, Δf_sinc + additional_offset; sample=sample_recovery)
        bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=sample_recovery)
        RRdelay = RR  - dur(bssfp) - dur(t2p) - dur(fatsat)
        seq += t2p
        seq += fatsat             
        seq += bssfp
        seq += sample_recovery ? ADC(40, RRdelay) : Delay(RRdelay) #Sampling recovery curve
    end
    return seq
end


sim_method = BlochDict(save_Mz=true)
#Sequence
seq = CMRA_iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, 
                        im_segments, iNAV_flip_angle, im_flip_angle; sample_recovery=true)
#Simulation
simParams = Dict{String,Any}("return_type"=>"mat", "sim_method"=>sim_method)
magnetization = simulate(obj, seq, sys; simParams)

# p, seqd = KomaMRIPlots.plot_seqd(seq; simParams); p
labels = ["Myocardium", "Blood", "Fat"]
colors = ["blue", "red", "green"]
spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)'] 
t = KomaMRICore.get_adc_sampling_times(seq) 

Mxy(i) = abs.(sum(magnetization[:,spin_group[i],1,1][:,1,:],dims=2)[:]/length(spin_group[i]))
Mz(i) = real.(sum(magnetization[:,spin_group[i],2,1][:,1,:],dims=2)[:]/length(spin_group[i]))
p0 = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
for i=eachindex(spin_group)
    p1 = scatter(x=t, y=Mxy(i), name=labels[i], legendgroup=labels[i], marker_color=colors[i])
    p2 = scatter(x=t, y=Mz(i), name=labels[i], legendgroup=labels[i], showlegend=false, marker_color=colors[i])#,line=attr(dash="dash"))
    add_trace!(p0, p1, row=1, col=1)
    add_trace!(p0, p2, row=2, col=1)
end
p0
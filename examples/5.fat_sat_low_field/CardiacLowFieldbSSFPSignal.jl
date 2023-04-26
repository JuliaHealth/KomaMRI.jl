using  KomaMRI, PlotlyJS, MAT
#Phantom, values from: 
# A review of normal tissue hydrogen NMR relaxation times and relaxation mechanisms from 1–100 MHz: 
# Dependence on tissue type, NMR frequency, temperature, species, excision, and age, Paul A. Bottomley,  
# Thomas H. Foster,  Raymond E. Argersinger,  Leah M. Pfeifer
B0 = 0.55
fat_ppm = 3.4e-6
Niso = 400
fat_freq = γ*B0*fat_ppm
off = range(-1, 1, Niso)
dx = 0.2*range(-1e-3, 1e-3, Niso)
myocard = Phantom{Float64}(x=dx, T1=834e-3*ones(Niso), T2=65e-3*ones(Niso))
blood =   Phantom{Float64}(x=dx, T1=960e-3*ones(Niso), T2=375e-3*ones(Niso))
fat =     Phantom{Float64}(x=dx, T1=92e-3*ones(Niso), T2=100e-3*ones(Niso), Δw=2π*(fat_freq .+ 5*off))
obj = myocard+blood+fat
#Sequence parameters
TR = 5.8e-3
TE = TR / 2 #bSSFP condition
Tadc = 1e-6
RR = 1 #1 [s]
Trf = 10e-6  #1 [ms]
B1 = 1 / (360*γ*Trf)
iNAV_lines = 6
im_segments = 20
iNAV_flip_angle = 3.2
im_flip_angle = 90
number_dummy_heart_beats = 3#10
#Pre-pulses
T2prep_duration = 50e-3
Tfatsat = 30e-3
Δf_sinc = fat_freq
additional_offset = 0
FatSat_flip_angle = 180
#Scanner
sys = Scanner()

function FatSat(α, Δf)
    RF_wf = matread("./examples/5.fat_sat_low_field/SINC_1.mat")["SINC_1"]
    seq = Sequence()
    seq += RF(RF_wf, Tfatsat, Δf)
    α_ref = get_flip_angles(seq)[2]
    seq *= (α/α_ref + 0im)  
    seq += ADC(1, 1e-6)
    seq += Grad(20e-3, 4.7e-3, 1e-3)
    seq += ADC(1, 1e-6)
    return seq
end

function T2prep(TE)
    seq = Sequence()
    seq += RF(90 * B1, Trf)
    seq += ADC(10, TE/2)
    seq += RF(180im * B1/2, 2Trf)
    seq += ADC(10, TE/2)
    seq += RF(-90 * B1, Trf)
    seq += ADC(1, 1e-6)
    seq += Grad(80e-3, 4.7e-3, 1e-3)
    seq += ADC(1, 1e-6)
    return seq
end

function bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
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
        seq += Delay(TE - Trf)
        seq += ADC(1, Tadc)
        seq += Delay(TR - TE - Tadc)
        k += 1
    end
    return seq
end

function iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
    k = 0
    seq = Sequence()
    for hb = 0 : number_dummy_heart_beats
        seq += T2prep(T2prep_duration)
        seq += FatSat(FatSat_flip_angle, Δf_sinc + additional_offset)
        seq += bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
        seq += ADC(40, RR - TR * (iNAV_lines + im_segments)) #Sampling recovery curve
    end
    return seq
end

seq_def = iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, 90); 
αs = [90.] #max.(1:10:180, iNAV_flip_angle)
sim_method = BlochDict(save_Mz=true)
magnetization = zeros(ComplexF64, length(αs), KomaMRICore.sim_output_dim(obj,seq_def,sys,sim_method)...) 
for (j, im_flip_angle) = enumerate(αs)
    #Sequence
    seq = iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
    #Simulation
    simParams = Dict{String,Any}("return_type"=>"mat", "sim_method"=>sim_method)
    magnetization[j,:,:,:] = simulate(obj, seq, sys; simParams)
end


labels = ["Myocardium", "Blood", "Fat"]
colors = ["blue", "red", "green"]
spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)'] 
t = KomaMRICore.get_adc_sampling_times(seq_def)
angle_idx = 1
p0 = plot_seq(seq_def;slider=false, show_sim_blocks=true)

Mxy(i) = abs.(sum(magnetization[angle_idx,:,spin_group[i],1],dims=3)[:,1]/length(spin_group[i]))
Mz(i) = real.(sum(magnetization[angle_idx,:,spin_group[i],2],dims=3)[:,1]/length(spin_group[i]))
p = make_subplots(rows=1, cols=1, shared_xaxes=true, vertical_spacing=0.02)
for i=eachindex(spin_group)
    p1 = scatter(x=t, y=Mxy(i),
            name=labels[i],legendgroup=labels[i],marker_color=colors[i])
    p2 = scatter(x=t, y=Mz(i),
            name=labels[i],legendgroup=labels[i], showlegend=false,marker_color=colors[i],line=attr(dash="dash"))
    add_trace!(p, p1, row=1, col=1)
    add_trace!(p, p2, row=1, col=1)
end
p
# kcenter_steady_state = number_dummy_heart_beats * RR + iNAV_lines * TR + TE
# p3 = plot(αs, abs.(sum(magnetization[:,argmin(abs.(t.-kcenter_steady_state))[1],:,1]), 
#     Layout(xaxis_title="Flip angle [deg]", title=attr(text= "Signal for different bSSFP at flip-angles")))
# if sim_method.save_Mz
#     display([p1; p2; p3])
# else
#     display([p1; p3])
# end
# display([p0; p])
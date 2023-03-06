using  KomaMRI, PlotlyJS
#Phantom, values from: 
# A review of normal tissue hydrogen NMR relaxation times and relaxation mechanisms from 1–100 MHz: 
# Dependence on tissue type, NMR frequency, temperature, species, excision, and age, Paul A. Bottomley,  
# Thomas H. Foster,  Raymond E. Argersinger,  Leah M. Pfeifer
B0 = 0.55
fat_ppm = 3.4e-6
Niso = 20
# myocard = Phantom{Float64}(x=[0], ρ=[1.0], T1=[600e-3], T2=[55e-3])  #50, 1333
# blood =   Phantom{Float64}(x=[0], ρ=[1.0], T1=[1200e-3],T2=[243e-3]) #243, 1489
# fat =     Phantom{Float64}(x=[0], ρ=[1.0], T1=[217e-3], T2=[95e-3], Δw=[-2π*γ*B0*fat_ppm])

myocard = Phantom{Float64}(x=[0], T1=[600e-3], T2=[60e-3])
blood =   Phantom{Float64}(x=[0], T1=[1200e-3],T2=[410e-3])
fat =     Phantom{Float64}(x=[0], T1=[217e-3], T2=[95e-3], Δw=[-2π*γ*B0*fat_ppm])
obj = myocard+blood+fat
#Sequence parameters
TR = 5.81e-3
TE = TR / 2 #bSSFP condition
Tadc = 1e-6
RR = 1 #1 [s]
Trf = 1e-4  #1 [ms]
B1 = 1 / (360*γ*Trf)
iNAV_lines = 14
im_segments = 20
iNAV_flip_angle = 3.2
im_flip_angle = 90
number_dummy_heart_beats = 4
#Pre-pulses
T2prep_duration = 60e-3
Tfatsat = 20e-3
FatSat_flip_angle = 180
#Scanner
sys = Scanner()

function FatSat(FA)
    sinc = PulseDesigner.RF_sinc(1 / (360*γ*Tfatsat), Tfatsat, sys; Δf=-γ*B0*fat_ppm)[1]
    α_ref = KomaMRICore.get_flip_angles(sinc)[1]
    sinc *= (FA/α_ref+0im)
end

function T2prep(TE)
    seq = Sequence()
    seq += RF(90 * B1, Trf)
    seq += Delay(TE/2-Trf)
    seq += RF(180im * B1, Trf)
    seq += Delay(TE/2-Trf)
    seq += RF(-90 * B1, Trf)
end

function iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
    k = 0
    seq = Sequence()
    for hb = 0 : number_dummy_heart_beats
        # seq += T2prep(T2prep_duration)
        # seq += FatSat(FatSat_flip_angle)
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
        # seq += Delay(RR - TR * (iNAV_lines + im_segments))
        seq += ADC(20, RR - TR * (iNAV_lines + im_segments)) #Sampling recovery curve
    end
    return seq
end

seq_def = iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, 90); plot_seq(seq_def)
αs = max.(1:10:180, iNAV_flip_angle)
sim_method = BlochDict(save_Mz=true)
magnetization = zeros(ComplexF64, length(αs), KomaMRICore.sim_output_dim(obj,seq_def,sys,sim_method)...) 
for (j, im_flip_angle) = enumerate(αs)
    #Sequence
    seq = iNAV_bSSFP_cardiac(number_dummy_heart_beats, iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
    #Simulation
    simParams = Dict{String,Any}("return_type"=>"mat", "sim_method"=>sim_method)
    magnetization[j,:,:,:] = simulate(obj, seq, sys; simParams)
end


labels= ["Myocardium", "Blood", "Fat"]
t = KomaMRICore.get_adc_sampling_times(seq_def)
# p0 = plot_seq(seq_def;slider=false)
p1 = plot([
    scatter(x=t, y=abs.(magnetization[findall(αs.==91)[1],:,i,1]),name=labels[i]) 
    for i=1:length(obj)], Layout(xaxis_title="Time [s]", title=attr(text= "Mxy evolution bSSFP at 0.55T with α = 91 deg"))
)
p2 = plot([
    scatter(x=t, y=abs.(magnetization[findall(αs.==91)[1],:,i,2]),name=labels[i]) 
    for i=1:length(obj)], Layout(xaxis_title="Time [s]", title=attr(text= "Mz evolution bSSFP at 0.55T with α = 91 deg"))
)
kcenter_steady_state = number_dummy_heart_beats * RR + iNAV_lines * TR + TE
p3 = plot(αs, abs.(magnetization[:,argmin(abs.(t.-kcenter_steady_state))[1],:,1]), Layout(xaxis_title="Flip angle [deg]", title=attr(text= "Signal for different bSSFP at flip-angles")))

if sim_method.save_Mz
    display([p1; p2; p3])
else
    display([p1; p3])
end
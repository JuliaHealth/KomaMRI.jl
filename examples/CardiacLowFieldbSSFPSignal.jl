using  KomaMRI
#Phantom
# B0 = 0.55
# fat_ppm = 3.4e-6
blood =   Phantom(x=[0],T1=[1200e-3],T2=[410e-3])
# fat =     Phantom(x=[0],T1=[1200e-3],T2=[410e-3],Δw=[2π*γ*B0*fat_ppm])
# myocard = Phantom(x=[0],T1=[1200e-3],T2=[410e-3])
#Sequence parameters
TR = 5.81e-3
TE = TR / 2 #bSSFP condition
RR = 1 #1 [s]
Trf = 1e-3  #1ms
Tadc = 3e-3 #1ms
B1 = 1 / (360*γ*Trf)
iNAV_lines = 14
im_segments = 19
first_flip_angle = 3.2 #3, 20
im_flip_angle = 90.0
number_heart_beats = 1
#Sequence loop
seq = Sequence()
for hb = 0:number_heart_beats - 1
    for i = 0:iNAV_lines + im_segments - 1
        α = min( (im_flip_angle - first_flip_angle) / iNAV_lines * i + first_flip_angle, im_flip_angle) * exp(1im * π * i)
        global seq += α * RF(Trf, B1)
        global seq += Delay(TE - Tadc / 2)
        global seq += ADC(10, Tadc)
        global seq += Delay(TR - TE - Tadc / 2)
    end
    global seq += Delay(RR - TR * (iNAV_lines + im_segments))
end
plot_seq(seq)

#Simulation
sys = Scanner()
raw = simulate(blood, seq, sys; simParams=Dict{String,Any}("Nblocks"=>1))
plot_signal(raw)
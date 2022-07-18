using KomaMRI, JLD2
#Scanner
sys = Scanner()
sys.ADC_Δt = 4e-6 #ADC sampling time
#Seq
α, T = π/2, 1e-3
B1 = α / (2π*γ*T) # α =  2π γ B1 T -> B1 = α / (2π γ T)
EX = PulseDesigner.RF_hard(B1,T,sys)
FOV, N = 25.6e-2, 101
EPI,_,_,_ = PulseDesigner.EPI_base(FOV,N,sys)
TE = 50e-3
d1 = TE - dur(EPI)/2 - dur(EX)
D1 = Delay(d1) #You can obtain the duration of the EPI with dur(EPI)
R = rotz(0.) #Play with the rotation of the acquisition!
# For a GE acquisition
seq = EX + D1 + R*EPI
#Save seq
println(seq)
@save "./myEPI_GE.seq" seq
# For a SE, the delays must refocuse at k=0 (half the EPI)
REF = 2.0im*EX #twice the B1, 90 deg -> 180 deg
halfTE = TE / 2
d1 = halfTE - dur(EX)
d2 = halfTE - dur(EPI)/2 - dur(REF)
D1 = Delay(d1)
D2 = Delay(d2)
seq = EX + D1 + REF + D2 + R*EPI
#Save seq
println(seq)
@save "./myEPI_SE.seq" seq
#Go to SpinLab GUI, loag "myEPI_GE.seq", press "Run"

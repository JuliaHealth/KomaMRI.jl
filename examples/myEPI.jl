using MRIsim, JLD2
using MRIsim: γ #<- [Hz/T] gamma barra
#Seq
α, T = π/2, 1e-3
B1 = α / (2π*γ*T) # α =  2π γ B1 T -> B1 = α / (2π γ T)
EX = PulseDesigner.RF_hard(B1,T)
FOV, N, Δt, Gmax = 25.6e-2, 101, 4e-6, 30e-3
EPI,_,_,_ = PulseDesigner.EPI_base(FOV,N,Δt,Gmax)
TE = 50e-3
d1 = TE - dur(EPI)/2 - dur(EX)
D1 = Sequence([delay(d1)]) #You can obtain the duration of the EPI with dur(EPI)
R = MRIsim.rotz(0.) #Play with the rotation of the acquisition!
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
D1 = Sequence([delay(d1)])
D2 = Sequence([delay(d2)])
seq = EX + D1 + REF + D2 + R*EPI
#Save seq
println(seq)
@save "./myEPI_SE.seq" seq
#Go to SpinLab GUI, loag "myEPI_GE.seq", press "Run"
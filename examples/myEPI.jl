using MRIsim, JLD2
using MRIsim: γ #<- [Hz/T] gamma barra
#Seq
α, T = π/2, 1e-3
B1 = α / (2π*γ*T) # α =  2π γ B1 T -> B1 = α / (2π γ T)
EX = PulseDesigner.RF_hard(B1,T)
FOV, N, Δt, Gmax = 25.6e-2, 101, 4e-6, 30e-3
EPI,_,_,_ = PulseDesigner.EPI_base(FOV,N,Δt,Gmax)
TE = 50e-3
d1 = Sequence([delay(TE - dur(EPI)/2 - dur(EX))]) #You can obtain the duration of the EPI with dur(EPI)
R = MRIsim.rotz(0.) #Play with the rotation of the acquisition!
# For a GE acquisition
seq = EX + d1 + R*EPI
#Save seq
println(seq)
@save "C:/Users/56988/Downloads/ProyectoIBM2101/myEPI_GE.seq" seq
# For a SE, the delays must refocuse at k=0 (half the EPI)
REF = 2.0im*EX #twice the B1, 90 deg -> 180 deg
halfTE = TE / 2
d1 = Sequence([delay(halfTE - dur(EX))])
d2 = Sequence([delay(halfTE - dur(EPI)/2 - dur(REF))])
seq = EX + d1 + REF + d2 + R*EPI
#Save seq
println(seq)
@save "C:/Users/56988/Downloads/ProyectoIBM2101/myEPI_SE.seq" seq
#Go to SpinLab GUI, loag "myEPI_GE.seq", press "Run"
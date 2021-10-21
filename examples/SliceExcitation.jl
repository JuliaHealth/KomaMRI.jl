# Example of rotation of acquisition and excitation
using MRIsim
using JLD2
# SEQ design
γ = 42.58e6
B1 = 6e-6; durRF = π/(2π*γ*B1)
EX = PulseDesigner.RF_hard(B1, durRF; G=[0 0])
Gmax = 60e-3
EPI,_,_,_ = PulseDesigner.EPI_base(40/100, 99, 4e-6, Gmax)
TE = 25e-3
d = delay(TE-dur(EPI)/2-dur(EX))
DELAY = Sequence([d;d])
#Rotation
Rz = MRIsim.rotz(-π/2) #counter-clockwise rotation
Ry = MRIsim.roty(-π/2) #counter-clockwise rotation
Rx = MRIsim.rotx(π/2)  #counter-clockwise rotation

seq = EX + DELAY + EPI*(Ry*Rz)
@save "./EPI_example.seq" seq=seq
plot_seq(seq)


# Example of rotation of acquisition and excitation
using KomaMRI
using JLD2
# Scanner
sys = Scanner()
sys.ADC_Δt = 4e-6
# SEQ design
B1 = 6e-6; durRF = π/(2π*γ*B1)
EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
EPI,_,_,_ = PulseDesigner.EPI_base(40/100, 101, sys)
TE = 25e-3
DELAY = Delay(TE-dur(EPI)/2-dur(EX))
#Rotation
Rz = KomaMRI.rotz(-π/2) #counter-clockwise rotation
Ry = KomaMRI.roty(-π/2) #counter-clockwise rotation
Rx = KomaMRI.rotx(π/2)  #counter-clockwise rotation

seq = EX + DELAY + EPI*(Ry*Rz)
@save "./EPI_example.seq" seq=seq
plot_seq(seq)


# Example of rotation of acquisition and excitation
using MRIsim
using JLD2
# SEQ design
γ = 42.58e6
B1 = 6e-6; durRF = π/2/(2π*γ*B1)
EX = PulseDesigner.RF_hard(B1, durRF; G=[2e-3 0])
Gmax = 60e-3
EPI,_,_,_ = PulseDesigner.EPI_base(40/100, 99, 4e-6, Gmax)
TE = 25e-3 
d = delay(TE-dur(EPI)/2-dur(EX))
DELAY = Sequence([d;d])
#Rotation
for t = 0:45:180
    Rz = MRIsim.rotz(t/180*π) #counter-clockwise rotation
    seq = EX*Rz + DELAY + EPI
    # plot_seq(seq)
    @save "./EPI_example_$t.seq" seq=seq
end


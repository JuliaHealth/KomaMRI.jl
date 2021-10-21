# Code to asses the effect of moment-compensation on the diffusion signal

## Read
#Waveform
using MRIsim, PlotlyJS
b = [0, 50, 200, 400, 600, 790]
E = Array{Any, 2}(undef,6,4)
FID = false
path_file = "/home/ccp/Documents/MRIsim.jl/"
#Timings
f = readlines(path_file*"qte_timings.txt")
timings = tryparse.(Float64, split(f[2]," "))
δ1, rf180, δ2 = timings[timings .≠ nothing]*1e-3
#Normal DTI
G, n1, n2 = MRIsim.read_diff_fwf(path_file*"QTE_Waveforms/MX_MC2_b790.txt")
Gmax = 21e-3
g1, g2 = G[1:20:n1,1], G[1:20:n2,4]
seq_dti = Sequence(Gmax*[Grad(g1,δ1) delay(rf180) -Grad(g2,δ2)])
#ACQ
if FID
    seq_ACQ = Sequence([Grad(0,0)], [RF(0,0)], DAC(1,0))
else
    FOV = 4e-2 #10cm
    Nxy = 20
    ΔT = 20e-6 #100um
    Gmax_acq = 60e-3
    seq_ACQ,_,_,_ = PulseDesigner.EPI_base(FOV, Nxy, ΔT, Gmax_acq)
end
#Phantom
Ns = 10_000; Nsx = floor(Int64, sqrt(Ns) ); Ns = Nsx*Nsx
ADC = 2e-9 #Normal diffusion coefficient
v = 6e-2 # 2mm/s Capilar flow
Δx = .5e-2 #3mm
x = range(-Δx/2,Δx/2,length=Nsx) #spin coordinates
y = range(-Δx/2,Δx/2,length=Nsx)  #spin coordinates
xx, yy = x .+ y'*0, x*0 .+ y' #grid points
T = dur(seq_dti)+dur(seq_ACQ)/2
phantom  = Phantom(x=xx[:], y=yy[:], Dλ1=ADC*ones(Ns), Dλ2=ADC*ones(Ns))
phantom2 = Phantom(x=xx[:], y=yy[:], ux=(x,y,z,t)->v*t, Dλ1=ADC*ones(Ns), Dλ2=ADC*ones(Ns))
#Iter
for idx = 1:5
#Seq
G, n1, n2 = MRIsim.read_diff_fwf(path_file*"QTE_Waveforms/MX_MC2_b$(b[idx+1]).txt")
g1, g2 = G[1:20:n1,1], G[1:20:n2,4]
seq = Sequence(Gmax*[Grad(g1,δ1) delay(rf180) -Grad(g2,δ2)])
seq += seq_ACQ
#Simulate
f = sqrt(b[idx+1]/790)
simParams = Dict(:Δt => 10e-6, :Nblocks=> 1)
if FID
    recParams = Dict(:skip_rec => true) 
else
    recParams = Dict(:epi => true, :recon => :fft, :Nx => Nxy)
end
S_MM = MRIsim.simulate(phantom, seq, simParams, recParams)
S_DTI = MRIsim.simulate(phantom, f*seq_dti+seq_ACQ, simParams, recParams)
S_MM2 = MRIsim.simulate(phantom2, seq, simParams, recParams)
S_DTI2 = MRIsim.simulate(phantom2, f*seq_dti+seq_ACQ, simParams, recParams)
S0 = MRIsim.simulate(phantom, 0*seq, simParams, recParams)
S02 = MRIsim.simulate(phantom2, 0*seq, simParams, recParams)
#Saving results
E[idx+1,1] = S_MM
E[idx+1,2] = S_DTI
E[idx+1,3] = S_MM2
E[idx+1,4] = S_DTI2
E[1,1] = S0
E[1,2] = S0
E[1,3] = S02
E[1,4] = S02
end

##Plot
i = 4
zmax = 12
plot(heatmap(z=abs.(E[i,1]), transpose = false, zmin=0, zmax=zmax), Layout(;title="No mov. MM-DTI b$(b[i])"))
plot(heatmap(z=abs.(E[i,2]), transpose = false, zmin=0, zmax=zmax), Layout(;title="No mov. DTI b$(b[i])"))
plot(heatmap(z=abs.(E[i,3]), transpose = false, zmin=0, zmax=zmax), Layout(;title="Mov. MM-DTI b$(b[i])"))
plot(heatmap(z=abs.(E[i,4]), transpose = false, zmin=0, zmax=zmax), Layout(;title="Mov. DTI b$(b[i])"))

p = MRIsim.plot_grads_moments(seq_dti+seq_ACQ)
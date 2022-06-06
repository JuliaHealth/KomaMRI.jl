# Code to asses the effect of moment-compensation on the diffusion signal
# MOST OF THE SYNTAX IN THIS FILE HAS CHANGED, RUN WITH CAUTION
## Read
#Waveform
using KomaMRI, PlotlyJS
b = [0, 50, 200, 400, 600, 790]
E = Array{Any, 2}(undef,6,4)
FID = false
path_file = "/home/ccp/Documents/KomaMRI.jl/"
#Timings
f = readlines(path_file*"qte_timings.txt")
timings = tryparse.(Float64, split(f[2]," "))
δ1, rf180, δ2 = timings[timings .≠ nothing]*1e-3
#Normal DTI
G, n1, n2 = KomaMRI.read_diff_fwf(path_file*"QTE_Waveforms/MX_MC2_b790.txt")
Gmax = 21e-3
g1, g2 = G[1:20:n1,1], G[1:20:n2,4]
seq_dti = Sequence(Gmax*[Grad(g1,δ1) Delay(rf180) -Grad(g2,δ2)])
#ACQ
if FID
    seq_ACQ = Sequence([Grad(0,0)], [RF(0,0)], ADC(1,0))
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
v = 1e-2 # 2mm/s Capilar flow
Δx = 1e-2 #3mm
x = range(-Δx/2,Δx/2,length=Nsx) #spin coordinates
y = range(-Δx/2,Δx/2,length=Nsx)  #spin coordinates
xx, yy = x .+ y'*0, x*0 .+ y' #grid points
T = dur(seq_dti)+dur(seq_ACQ)/2
phantom  = Phantom(x=xx[:], y=yy[:], Dλ1=ADC*ones(Ns), Dλ2=ADC*ones(Ns))
phantom2 = Phantom(x=xx[:], y=yy[:], ux=(x,y,z,t)-> (x .> Δx*.4).*(-v*t) .+ (x .< -Δx*.4).*(v*t), Dλ1=ADC*ones(Ns), Dλ2=ADC*ones(Ns))
#Iter
for idx = [2]
#Seq
G, n1, n2 = KomaMRI.read_diff_fwf(path_file*"QTE_Waveforms/MX_MC2_b$(b[idx+1]).txt")
g1, g2 = G[1:20:n1,1], G[1:20:n2,4]
seq = Sequence(Gmax*[Grad(g1,δ1) Delay(rf180) -Grad(g2,δ2)])
seq += seq_ACQ
#Simulate
f = sqrt(b[idx+1]/790)
simParams = Dict(:Δt => 10e-6, :Nblocks=> 1)
if FID
    recParams = Dict(:skip_rec => true) 
else
    recParams = Dict(:epi => true, :recon => :fft, :Nx => Nxy)
end
S_MM = simulate(phantom, seq; simParams) #Eventualmente se agregará el objeto hardware
S_DTI = simulate(phantom, f*seq_dti+seq_ACQ; simParams)
S_MM2 = simulate(phantom2, seq; simParams)
S_DTI2 = simulate(phantom2, f*seq_dti+seq_ACQ; simParams)
S0 = simulate(phantom, 0*seq; simParams)
S02 = simulate(phantom2, 0*seq; simParams)
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

## Plot
i = 3
zmax = 15
plot(heatmap(z=abs.(E[i,1]), transpose=false, zmin=0, zmax=zmax), Layout(;title="No mov. MM-DTI b$(b[i])"))
plot(heatmap(z=abs.(E[i,2]), transpose=false, zmin=0, zmax=zmax), Layout(;title="No mov. DTI b$(b[i])"))
plot(heatmap(z=abs.(E[i,3]), transpose=false, zmin=0, zmax=zmax), Layout(;title="Mov. MM-DTI b$(b[i])"))
plot(heatmap(z=abs.(E[i,4]), transpose=false, zmin=0, zmax=zmax), Layout(;title="Mov. DTI b$(b[i])"))

p = KomaMRI.plot_grads_moments(seq_dti+seq_ACQ)
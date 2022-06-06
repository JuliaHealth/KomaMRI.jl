using Koma, MAT, JLD2
using MRIReco
#Phantom 
phantom = brain_phantom2D()
phantom.Δw .*= 0 #Removing off-resonance
#Scanner
sys = Scanner()
#Seq: Excitation and radial acquisition
TB1 = 1e-3
EX = PulseDesigner.RF_hard(1/(2π*γ*TB1),TB1,sys) #α = γ B1 T
FOV, Nr = 23e-2, 100
RAD = PulseDesigner.radial_base(FOV, Nr, sys)
#MRF parameters
path = @__DIR__
data = matread(path*"/MRFinput.mat")
α, ϕ, TE, TR = data["FA"]/180*π, data["PH"]/180*π, data["TE"]*1e-3, data["TR"]*1e-3
A = α.*exp.(1im*ϕ)
#Inversion pulse
TINV = 50e-3
INV = (π+0im)*EX + Delay(TINV)
#Delays
Ta = dur(RAD)
delayTE = [Delay(TE[n]-Ta) for n=1:length(α)]
delayTR = [Delay(TR[n]-TE[n]-Ta) for n=1:length(α)]
#MRF with rotated spokes
φ, Nφ = (√5 + 1)/2, 7; 
Δθ = π/(φ+Nφ-1) # Tiny golden angle with Nφ = 7
NTRs = 200 #Number of TRs
seq = INV + sum([A[n]*EX + delayTE[n] + rotz((n-1)*Δθ)*RAD + delayTR[n] for n=1:NTRs])
seq.DEF["Nz"] = NTRs
jldsave("./examples/3.mrf_ismrm/mrf.seqk"; seq=seq)
# plot_seq(seq)
## Simulation
fingerprint = simulate(phantom, seq, sys); #This takes like 5 min for NTRs = 500.
## Output ISMRMRD
raw_ismrmrd = KomaMRI.rawSignalToISMRMRD([fingerprint;;],seq;phantom,sys)
fname = "MRF_signal_Δθ_$(floor(Int64, Δθ/π*180))_NTRs_$NTRs"
fout = ISMRMRDFile("./examples/3.mrf_ismrm/$fname.h5")
save(fout, raw_ismrmrd)
## Output MAT
# matwrite("./examples/3.mrf_ismrm/$fname.mat", Dict("signal"=>fingerprint))
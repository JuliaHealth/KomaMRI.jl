using MRIsim, MAT, PlotlyJS
using MRIsim: γ, rotz, brain_phantom2D
#Phantom 
phantom = brain_phantom2D()
phantom.Δw = zeros(size(phantom)) #Removing off-resonance
#Seq: Excitation and radial acquisition
TB1 = 1e-3
EX = PulseDesigner.RF_hard(1/(2π*γ*TB1),TB1) #α = γ B1 T
FOV, Nr, T, Gmax = 46e-2, 100, 50e-6, 30e-3
RAD, _, Δθ, Ta = PulseDesigner.radial_base(FOV, Nr, T, Gmax)
#MRF parameters
data = matread("./src/data/MRFinput.mat")
α, ϕ, TE, TR = data["FA"]/180*π, data["PH"]/180*π, data["TE"]*1e-3, data["TR"]*1e-3
A = α.*exp.(1im*ϕ)
#Inversion pulse
TINV = 50e-3
INV = (π+0im)*EX + Sequence([delay(TINV)])
#Delays
delayTE = [Sequence([delay(TE[n]-Ta)]) for n=1:length(α)]
delayTR = [Sequence([delay(TR[n]-TE[n]-Ta)]) for n=1:length(α)]
#MRF with rotated spokes
# φ, Nφ = (√5 + 1)/2, 7; Δθ = π/(φ+Nφ-1) # Uncomment for tiny golden angle 7
NTRs = 1000 #Number of TRs
seq = INV + sum([A[n]*EX + delayTE[n] + rotz((n-1)*Δθ)*RAD + delayTR[n] for n=1:NTRs])
plot_seq(seq)
## Simulation
simParams = Dict(:step => "uniform", :Δt => 25e-6, :Nblocks=> 20*NTRs)
recParams = Dict(:skip_rec => true)
fingerprint = MRIsim.simulate(phantom, seq, simParams, recParams)
plot(abs.(fingerprint))
matwrite("./MRF_signal_Δθ_$(floor(Int64, Δθ/π*180)).mat", Dict("signal"=>fingerprint))
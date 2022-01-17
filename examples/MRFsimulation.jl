using Koma, MAT, PlotlyJS
#Phantom 
phantom = brain_phantom2D()
phantom.Δw = zeros(size(phantom)) #Removing off-resonance
#Scanner
sys = Scanner()
#Seq: Excitation and radial acquisition
TB1 = 1e-3
EX = PulseDesigner.RF_hard(1/(2π*γ*TB1),TB1) #α = γ B1 T
FOV, Nr = 46e-2, 100
RAD = PulseDesigner.radial_base(FOV, Nr, sys)
#MRF parameters
data = matread("./src/examples/MRFinput.mat")
α, ϕ, TE, TR = data["FA"]/180*π, data["PH"]/180*π, data["TE"]*1e-3, data["TR"]*1e-3
A = α.*exp.(1im*ϕ)
#Inversion pulse
TINV = 50e-3
INV = (π+0im)*EX + Sequence([Delay(TINV)])
#Delays
delayTE = [Delay(TE[n]-Ta) for n=1:length(α)]
delayTR = [Delay(TR[n]-TE[n]-Ta) for n=1:length(α)]
#MRF with rotated spokes
# φ, Nφ = (√5 + 1)/2, 7; Δθ = π/(φ+Nφ-1) # Uncomment for tiny golden angle 7
NTRs = 1000 #Number of TRs
seq = INV + sum([A[n]*EX + delayTE[n] + rotz((n-1)*Δθ)*RAD + delayTR[n] for n=1:NTRs])
plot_seq(seq)
## Simulation
simParams = Dict("step" => "uniform", "Δt" => 25e-6, "Nblocks"=> 20*NTRs)
fingerprint = simulate(phantom, seq, scanner, simParams)
plot(abs.(fingerprint))
matwrite("./MRF_signal_Δθ_$(floor(Int64, Δθ/π*180)).mat", Dict("signal"=>fingerprint))
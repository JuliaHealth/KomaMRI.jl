# Compare speed with MRiLab
using KomaMRI
#Scanner object
sys = Scanner()
sys.Smax=150    # [mT/m/ms]
sys.Gmax=500e-3 # [T/m]
FOV = 0.2       # [m]
N = 80          # Reconstructed image N×N
#RF sinc
B1 = 24.7835e-6 # For 90 deg flip angle
Trf = 1e-3
rf = PulseDesigner.RF_sinc(B1, Trf, sys; TBP=4)
#Spiral sequence
TE = 50e-3  # 50e-3 [s]
TR = 10     # 10 [s]
Nint = 8
spiral = PulseDesigner.spiral_base(FOV, N, sys; λ=2.1, BW=120e3, Nint)
delayTE = Delay(TE-Trf/2)
delayTR = Delay(TR-Trf/2-TE-dur(spiral(0)))
seq = Sequence()
for i = 1:Nint
    global seq += rf + delayTE + spiral(i-1) + delayTR
end
plot_seq(seq; range=[0 TE+dur(spiral(0))+1e-3]*1e3)
plot_kspace(seq)
# Phantom
filepath = "examples/3.koma_paper/comparison_mrilab/phantom/"
filename = filepath * "brain_mrilab.mat"
FRange_filename = filepath * "FRange.mat"
phantom = KomaMRI.read_phantom_MRiLab(filename; FRange_filename)
# Simulate
simParams=Dict{String,Any}("Nblocks"=>1, "gpu"=>true, "Nthreads"=>1, "return_type"=>"raw", "Δt_rf"=>5e-5)
raw = simulate(phantom, seq, sys; simParams)
using CUDA
CUDA.@profile simulate(phantom, seq, sys; simParams);
plot_signal(raw) #; show_sim_blocks=true)
# # Recon
acq = AcquisitionData(raw)
reconParams = Dict{Symbol,Any}(
    :reco=>"direct",
    :reconSize=>(N, N),
    :densityWeighting => true
    )
image = reconstruction(acq, reconParams)
plot_image(abs.(image[:, :, 1]); height=400, width=400)
# Compare speed with MRiLab
using KomaMRI
## Scanner
sys = Scanner()
sys.Smax = 150    # [mT/m/ms]
sys.Gmax = 500e-3 # [T/m]
sys.GR_Δt = 4e-6  # [s]
FOV = 0.2       # [m]
N = 80          # Reconstructed image N×N
## Pulse programming
# RF sinc
B1 = 24.7835e-6 # For 90 deg flip angle
Trf = 1e-3
rf = PulseDesigner.RF_sinc(B1, Trf, sys; TBP=4)
# Spiral sequence
TE = 50e-3  # 50e-3 [s]
TR = 10     # 10 [s]
Nint = 8
spiral = PulseDesigner.spiral_base(FOV, N, sys; λ=2.1, BW=120e3, Nint)
delayTE = Delay(TE - Trf / 2)
delayTR = Delay(TR - Trf / 2 - TE - dur(spiral(0)))
seq = Sequence()
for i = 1:Nint
    global seq += rf + delayTE + spiral(i - 1) + delayTR
end
plot_seq(seq; range=[0 TE + dur(spiral(0)) + 1e-3] * 1e3)
plot_kspace(seq)
## Phantom
filepath = "examples/3.koma_paper/comparison_mrilab/phantom/"
filename = filepath * "brain_mrilab.mat"
FRange_filename = filepath * "FRange.mat"
phantom = read_phantom_MRiLab(filename; FRange_filename)
## Simulation
simParams = Dict{String,Any}(
    "Nblocks" => 20,
    "gpu" => true,
    "Nthreads" => 1
)

raw = simulate(phantom, seq, sys; simParams)
plot_signal(raw; range=[50.5, 54]) #; show_sim_blocks=true)

using CUDA
CUDA.@profile ( simulate(phantom, seq, sys; simParams) );
## Recon
acq = AcquisitionData(raw)
reconParams = Dict{Symbol,Any}(
    :reco => "direct",
    :reconSize => (N, N),
    :densityWeighting => true
)
image = reconstruction(acq, reconParams)
plot_image(abs.(image[:, :, 1]); height=380, width=400)
# To profile this code use:
# nsys launch julia --project=. examples/3.koma_paper/comparison_mrilab/MRiLab_speed.jl
# Then, open the report by using NVIDIA Nsight Systems with nsys-ui
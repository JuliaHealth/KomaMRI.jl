# Compare speed of multi-shot spiral with MRiLab
cd(dirname(@__FILE__))
using KomaMRI, Suppressor
## Scanner
sys = Scanner()
sys.Smax = 150    # [mT/m/ms]
sys.Gmax = 500e-3 # [T/m]
sys.GR_Δt = 4e-6  # [s]
FOV = 0.2         # [m]
N = 80            # Reconstructed image N×N
## Pulse programming
# RF sinc
B1 = 1e-6 # For 90 deg flip angle
Trf = 1e-3
rf = PulseDesigner.RF_sinc(B1, Trf, sys; TBP=4, a=0.46)
α_desired = 90 + 0im
α =  KomaMRI.get_flip_angles(rf)[1]
rf *= α_desired / α #Linearly adjusts B1 to obtain desired FA
# Spiral sequence
TE = 50e-3  # 50e-3 [s]
TR = 10     # 10 [s]
Nint = 8
spiral = PulseDesigner.spiral_base(FOV, N, sys; BW=60e3, Nint, λ=2.1)
delayTE = Delay(TE - Trf / 2)
delayTR = Delay(TR - Trf / 2 - TE - dur(spiral(0)))
seq = Sequence()
for i = 1:Nint
    global seq += rf + delayTE + spiral(i - 1) + delayTR
end
plot_seq(seq; range=[0 TE + dur(spiral(0)) + 1e-3] * 1e3)
plot_kspace(seq)
## Phantom
filepath = "./phantom/"
filename = filepath * "brain_mrilab.mat"
FRange_filename = filepath * "FRange.mat" #Slab within slice thickness
phantom = read_phantom_MRiLab(filename; FRange_filename)
## Simulation
if (ARGS == String[]) #No arguments, use defaults
    sim_params = Dict{String,Any}(
        "gpu" => true,
        "gpu_device" => 0
    )
else
    sim_params = Dict{String,Any}(
        "gpu" => ARGS[1] == "gpu" ? true : false,
        "gpu_device" => parse(Int64, ARGS[2])
    )
end

Nexp = 20
raw = @suppress simulate(phantom, seq, sys; sim_params) #warmup
for i = 1:Nexp
    local raw = simulate(phantom, seq, sys; sim_params)
end
# plot_signal(raw; range=[50.5, 54]) #; show_sim_blocks=true)

# To profile this code use:
# nsys launch julia --project=. --color=yes examples/5.koma_paper/comparison_mrilab/MRiLab_speed.jl
# Then, open the report by using NVIDIA Nsight Systems with nsys-ui
# using CUDA
# CUDA.@profile ( simulate(phantom, seq, sys; sim_params) );

## Recon
# acq = AcquisitionData(raw)
# reconParams = Dict{Symbol,Any}(
#     :reco => "direct",
#     :reconSize => (N, N),
#     :densityWeighting => true
# )
# image = reconstruction(acq, reconParams)
# plot_image(abs.(image[:, :, 1]); height=380, width=400)

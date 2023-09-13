# Import packages
using PlotlyJS, KomaMRI

# Define scanner and object
sys = Scanner()
obj = brain_phantom2D()

# Define the sequence
function create_seq_epi(sys)
    B1 = sys.B1;
    durRF = π/2/(2π*γ*B1)
    EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
    N = 101
    FOV = 23e-2
    EPI = PulseDesigner.EPI(FOV, N, sys)
    TE = 30e-3
    d1 = TE-dur(EPI)/2-dur(EX)
    if d1 > 0 DELAY = Delay(d1) end
    seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI
    seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
    return seq
end
seq = create_seq_epi(sys)

# Plot the sequence
function plotsequence(seq::Sequence)
    sq = sequencevalues(seq)
    prfa = scatter(;x=sq.t[sq.rf_onmask], y=abs.(sq.rfa)[sq.rf_onmask], mode="lines+markers", name="RF")
    pgxa = scatter(;x=sq.t[sq.gx_onmask], y=sq.gxa[sq.gx_onmask], mode="lines+markers", name="GX")
    pgya = scatter(;x=sq.t[sq.gy_onmask], y=sq.gya[sq.gy_onmask], mode="lines+markers", name="GY")
    pgza = scatter(;x=sq.t[sq.gz_onmask], y=sq.gza[sq.gz_onmask], mode="lines+markers", name="GZ")
    padc = scatter(;x=sq.tadc, y=(.01 .* ones(length(sq.tadc))), mode="markers", name="ADC")
    display(plot([prfa; pgxa; pgya; pgza; padc], Layout(title="New Discretized Sequence")))
end
plotsequence(seq)

# Plot the kspace of the sequence
function plotkspace(seq::Sequence)
    kx, ky, kz, adc_onmask = KomaMRICore.kspace(seq)
    p1 = scatter3d(; x=kx, y=ky, z=kz, mode="lines", hoverinfo="skip")
    p2 = scatter3d(; x=kx[adc_onmask], y=ky[adc_onmask], z=kz[adc_onmask], mode="markers")
    p3 = scatter3d(; x=[0.], y=[0.], z=[0.], marker=attr(symbol="cross",size=10,color="red"))
    display(plot([p1; p2; p3], Layout(title="K-Space of the New Discretized Sequence")))
end
plotkspace(seq)

# Simulate for new simple simulation
sq = sequencevalues(seq)
t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, tadc, blk_range = sq.t, sq.Δt, sq.rfa, sq.rfΔf, sq.gxa, sq.gya, sq.gza, sq.rf_onmask, sq.gx_onmask, sq.gy_onmask, sq.gz_onmask, sq.adc_onmask, sq.tadc, sq.blk_range
magxy, magz, sig = komasim(seq, obj)
display(plot([scatter(;x=t, y=abs.([sum(magxy; dims=1)...]), mode="lines+markers", name="sum(magxy)"); scatter(;x=t, y=abs.(sig), mode="lines+markers", name="sig")], Layout(title="Raw-Signal of the New Simulator-Function")))

# Simulate for old simulation
simParams = KomaMRICore.default_sim_params()
simParams["return_type"] = "mat"
simParams["gpu"] = false
told = KomaMRICore.get_adc_sampling_times(seq)
sigold = simulate(obj, seq, sys; simParams)
display(plot([scatter(;x=told, y=abs.(sigold[:,1,1]), mode="lines+markers", name="old")], Layout(title="Raw-Signal of the Old Simulator-Function")))

# Compare new and old simulations
display(plot([scatter(;x=t[adc_onmask], y=abs.(sig[adc_onmask]), mode="lines+markers", name="new"); scatter(;x=told, y=abs.(sigold[:,1,1]), mode="lines+markers", name="old")], Layout(title="Comparison of the Raw-Signals")))

# Reconstruction
function recon(sig, seq; isnew=true)
    sigm = isnew ? reshape(sig, length(sig), 1, 1) : sig
    raw = signal_to_raw_data(sigm, seq)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqData.traj[1].nodes[:])) #Normalize k-space to -.5 to .5 for NUFFT
    Nx, Ny = raw.params["reconSize"][1:2]
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    image3d  = reshape(rec.data, Nx, Ny, :)
    image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
    return image2d
end
image2d = recon(sig[adc_onmask], seq; isnew=true)
display(plot_image(image2d, zmin=minimum(image2d), zmax=maximum(image2d), title="Reconstruction for New Simulator-Function"))
image2dold = recon(sigold, seq; isnew=false)
display(plot_image(image2dold, zmin=minimum(image2dold), zmax=maximum(image2dold), title="Reconstruction for Old Simulator-Function"))


############################################################################################
### Comparison with JEMRIS #################################################################
############################################################################################
using PlotlyJS, KomaMRI, HDF5

# Define input objects
seq = read_seq(joinpath(dirname(dirname(pathof(KomaMRI))), "KomaMRICore", "test", "test_files", "epi_100x100_TE100_FOV230.seq"))
obj = read_phantom_jemris(joinpath(dirname(dirname(pathof(KomaMRI))), "KomaMRICore", "test", "test_files", "sphere_chemical_shift.h5"))
sys = Scanner()

# Plot the sequence
plot_seq(seq)
plotsequence(seq)

# Old Koma simulation
simParams = Dict{String, Any}("gpu"=>false, "Nthreads"=>1, "sim_method"=>KomaMRICore.Bloch(), "return_type"=>"mat")
sigold = (simulate(obj, seq, sys; simParams) / prod(size(obj)))[:,1,1]
told = KomaMRICore.get_adc_sampling_times(seq)

# New Koma simulation
sq = sequencevalues(seq)
t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, tadc, blk_range = sq.t, sq.Δt, sq.rfa, sq.rfΔf, sq.gxa, sq.gya, sq.gza, sq.rf_onmask, sq.gx_onmask, sq.gy_onmask, sq.gz_onmask, sq.adc_onmask, sq.tadc, sq.blk_range
magxy, magz, sig = komasim(seq, obj)
signew = sig[adc_onmask] / prod(size(obj))
tnew = t[adc_onmask]

# Jemris ground-true data
sigjem = h5open(joinpath(dirname(dirname(pathof(KomaMRI))), "KomaMRICore", "test", "test_files", "jemris_signals_epi_sphere_cs.h5"))["/signal/channels/00"]
sigjem = sigjem[1,:] + 1im*sigjem[2,:]
sigjem = sigjem[:]

# Compare signals in plots
display(plot([scatter(; y=abs.(sigold), mode="lines+markers", name="old"); scatter(; y=abs.(signew), mode="lines+markers", name="new");], Layout(title="Comparison of the Raw-Signals")))
display(plot([scatter(; y=abs.(sigjem), mode="lines+markers", name="jem"); scatter(; y=abs.(sigold), mode="lines+markers", name="old");], Layout(title="Comparison of the Raw-Signals")))
display(plot([scatter(; y=abs.(sigjem), mode="lines+markers", name="jem"); scatter(; y=abs.(signew), mode="lines+markers", name="new");], Layout(title="Comparison of the Raw-Signals")))

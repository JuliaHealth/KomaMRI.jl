# Import packages
using PlotlyJS, KomaMRI

# For defining an EPI sequence
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

# Define scanner, object and sequence
sys = Scanner()
seq = create_seq_epi(sys)
obj = brain_phantom2D()

# Define simulator parameters
Δtgr, Δtrf = 1e-3, 1e-5
simParams = KomaMRICore.default_sim_params()
simParams["Δt"], simParams["Δt_rf"] = Δtgr, Δtrf
simParams["Nthreads"], simParams["Nblocks"] = 1, 1
simParams["return_type"] = "mat"
simParams["gpu"] = false
simParams["precision"] == "f64"

# Get the discretized sequences
seqori = @time KomaMRICore.discretize(seq; simParams);
sqs = samples(seq, Δtgr, Δtrf)
seqnew = KomaMRICore.SEQD(sqs.Δt, sqs.t, complex.(sqs.rfa), sqs.rfΔfc, sqs.gxa, sqs.gya, sqs.gza, sqs.adconmask)

# Plot the discretized sequences
prfaori = scatter(; x=seqori.t, y=(5000 .* abs.(seqori.B1)), mode="lines+markers", name="RF")
pgxaori = scatter(; x=seqori.t, y=seqori.Gx, mode="lines+markers", name="GX")
pgyaori = scatter(; x=seqori.t, y=seqori.Gy, mode="lines+markers", name="GY")
pgzaori = scatter(; x=seqori.t, y=seqori.Gz, mode="lines+markers", name="GZ")
padcori = scatter(; x=seqori.t[seqori.ADC], y=(.01 .* ones(length(seqori.ADC))), mode="markers", name="ADC")
prfanew = scatter(; x=seqnew.t, y=(5000 .* abs.(seqnew.rfa)), mode="lines+markers", name="RF")
pgxanew = scatter(; x=seqnew.t, y=seqnew.gxa, mode="lines+markers", name="GX")
pgyanew = scatter(; x=seqnew.t, y=seqnew.gya, mode="lines+markers", name="GY")
pgzanew = scatter(; x=seqnew.t, y=seqnew.gza, mode="lines+markers", name="GZ")
padcnew = scatter(; x=seqnew.t[seqnew.adconmask], y=(.01 .* ones(length(seqnew.adconmask))), mode="markers", name="ADC")
display(plot([prfaori; pgxaori; pgyaori; pgzaori; padcori; prfanew; pgxanew; pgyanew; pgzanew; padcnew], Layout(title="Seq comparison")))

# Simulate for original and new discretization-and-simulator
tori = KomaMRICore.get_adc_sampling_times(seq)
sigori = @time simulate(obj, seq, sys; simParams)[:];   # this asumes that the raw signal have only one column component with data
signew = @time seqsim(seq, obj, Δtgr, Δtrf; onadc=true);
signewall = @time seqsim(seq, obj, Δtgr, Δtrf; onadc=false);

display(plot([scatter(; x=seqnew.t[seqnew.adconmask], y=abs.(signew), mode="lines+markers", name="test-sim"); scatter(;x=tori, y=abs.(sigori), mode="lines+markers", name="master-sim")], Layout(title="Comparison of the Raw-Signals")))
display(plot([scatter(; x=seqnew.t, y=abs.(signewall), mode="lines+markers", name="test-sim-all"); scatter(;x=tori, y=abs.(sigori), mode="lines+markers", name="master-sim")], Layout(title="Comparison of the Raw-Signals")))

# Reconstruction
function recon(sig, seq)
    sigm = sig
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

# Reconstruct for original discretization
image2dori = recon(sigori, seq)
display(plot_image(image2dori, zmin=minimum(image2dori), zmax=maximum(image2dori), title="Reconstruction for Master Simulator"))

# Reconstruct for test simulator and new discretization
image2dnew = recon(signew, seq)
display(plot_image(image2dnew, zmin=minimum(image2dnew), zmax=maximum(image2dnew), title="Reconstruction for Test Simulator"))


############################################################################################
### TEST COMPARISON WITH MATHEMATICA #######################################################
############################################################################################
Tadc = 1e-3
Trf = Tadc
T1 = 1000e-3
T2 = 20e-3
Δw = 2π * 100
B1 = 2e-6 * (Tadc / Trf)
N = 6
sys = Scanner()
obj = Phantom{Float64}(x=[0],T1=[T1],T2=[T2],Δw=[Δw])
seq = Sequence()
seq += ADC(N, Tadc)
for i=1:3
    seq += RF(B1, Trf)
    seq += ADC(N, Tadc)
end
#seq = seq[2:end]

# Get the discretized sequences
seqori = @time KomaMRICore.discretize(seq; simParams);
sqs = samples(seq, Δtgr, Δtrf)
seqnew = KomaMRICore.SEQD(sqs.Δt, sqs.t, complex.(sqs.rfa), sqs.rfΔfc, sqs.gxa, sqs.gya, sqs.gza, sqs.adconmask)

# Plot the discretized sequences
prfaori = scatter(; x=seqori.t, y=(5000 .* abs.(seqori.B1)), mode="lines+markers", name="RF")
pgxaori = scatter(; x=seqori.t, y=seqori.Gx, mode="lines+markers", name="GX")
pgyaori = scatter(; x=seqori.t, y=seqori.Gy, mode="lines+markers", name="GY")
pgzaori = scatter(; x=seqori.t, y=seqori.Gz, mode="lines+markers", name="GZ")
padcori = scatter(; x=seqori.t[seqori.ADC], y=(.01 .* ones(length(seqori.ADC))), mode="markers", name="ADC")
prfanew = scatter(; x=seqnew.t, y=(5000 .* abs.(seqnew.rfa)), mode="lines+markers", name="RF")
pgxanew = scatter(; x=seqnew.t, y=seqnew.gxa, mode="lines+markers", name="GX")
pgyanew = scatter(; x=seqnew.t, y=seqnew.gya, mode="lines+markers", name="GY")
pgzanew = scatter(; x=seqnew.t, y=seqnew.gza, mode="lines+markers", name="GZ")
padcnew = scatter(; x=seqnew.t[seqnew.adconmask], y=(.01 .* ones(length(seqnew.adconmask))), mode="markers", name="ADC")
display(plot([prfaori; pgxaori; pgyaori; pgzaori; padcori; prfanew; pgxanew; pgyanew; pgzanew; padcnew], Layout(title="Seq comparison")))

# Simulate for original and new discretization-and-simulator
tori = KomaMRICore.get_adc_sampling_times(seq)
sigori = @time simulate(obj, seq, sys; simParams)[:];
signew = @time seqsim(seq, obj, Δtgr, Δtrf; onadc=true);
signewall = @time seqsim(seq, obj, Δtgr, Δtrf; onadc=false);

# Compare simulations
display(plot([scatter(; x=seqnew.t[seqnew.adconmask], y=abs.(signew), mode="lines+markers", name="test-sim"); scatter(;x=tori, y=abs.(sigori), mode="lines+markers", name="master-sim")], Layout(title="Comparison of the Raw-Signals")))
display(plot([scatter(; x=seqnew.t, y=abs.(signewall), mode="lines+markers", name="test-sim-all"); scatter(;x=tori, y=abs.(sigori), mode="lines+markers", name="master-sim")], Layout(title="Comparison of the Raw-Signals")))


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

# Old Koma simulation
sigold = (simulate(obj, seq, sys; simParams) / prod(size(obj)))[:,1,1]
told = KomaMRICore.get_adc_sampling_times(seq)

# New Koma simulation
sq = samples(seq, Δtgr, Δtrf)
sig = seqsim(seq, obj, Δtgr, Δtrf)
signew = sig[sq.adconmask] / prod(size(obj))
tnew = sq.t[sq.adconmask]

# Jemris ground-true data
sigjem = h5open(joinpath(dirname(dirname(pathof(KomaMRI))), "KomaMRICore", "test", "test_files", "jemris_signals_epi_sphere_cs.h5"))["/signal/channels/00"]
sigjem = sigjem[1,:] + 1im*sigjem[2,:]
sigjem = sigjem[:]

# Compare signals in plots
display(plot([scatter(; y=abs.(sigold), mode="lines+markers", name="old"); scatter(; y=abs.(signew), mode="lines+markers", name="new");], Layout(title="Comparison of the Raw-Signals")))
display(plot([scatter(; y=abs.(sigjem), mode="lines+markers", name="jem"); scatter(; y=abs.(sigold), mode="lines+markers", name="old");], Layout(title="Comparison of the Raw-Signals")))
display(plot([scatter(; y=abs.(sigjem), mode="lines+markers", name="jem"); scatter(; y=abs.(signew), mode="lines+markers", name="new");], Layout(title="Comparison of the Raw-Signals")))

#grs = grsamples(seq, 1; addblklim=false, addseqfirst=false, addseqlast=false, addblkfirst=false, addblklast=false, addrfx=false, addadc=false)
grs = grsamples(seq; addseqfirst=true, addseqlast=true, addblkfirst=false, addblklast=false, addrfx=true, addadc=true)
display(plot([scatter(; x=grs.t, y=grs.ax, mode="lines+markers", name="GX")]))

rfs = rfsamples(seq; addblklim=true, addseqfirst=false, addseqlast=false, addblkfirst=false, addblklast=false, addrfx=false)
display(plot([scatter(; x=rfs.t, y=rfs.a, mode="lines+markers", name="RF")]))

k = KomaMRICore.kspace(seq)
p1 = PlotlyJS.scatter3d(x=k.x, y=k.y, z=k.z, mode="lines", name="Trajectory", hoverinfo="skip")
p2 = PlotlyJS.scatter3d(x=k.x[k.adconmask], y=k.y[k.adconmask], z=k.z[k.adconmask], mode="markers", marker=attr(size=2), name="ADC")
p3 = PlotlyJS.scatter3d(x=[0], y=[0], z=[0], name="k=0", marker=attr(symbol="cross",size=10,color="red"))
plot([p1, p2, p3])

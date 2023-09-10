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
t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, adct = sequencevalues(seq)
prfa = scatter(;x=t[rf_onmask], y=abs.(rfa)[rf_onmask], mode="lines+markers", name="RF")
pgxa = scatter(;x=t[gx_onmask], y=gxa[gx_onmask], mode="lines+markers", name="GX")
pgya = scatter(;x=t[gy_onmask], y=gya[gy_onmask], mode="lines+markers", name="GY")
pgza = scatter(;x=t[gz_onmask], y=gza[gz_onmask], mode="lines+markers", name="GZ")
padc = scatter(;x=adct, y=(.01 .* ones(length(adct))), mode="markers", name="ADC")
display(plot([prfa; pgxa; pgya; pgza; padc]))

# Plot the kspace of the sequence
kx, ky, kz, adc_onmask = KomaMRICore.kspace(seq)
p1 = scatter3d(; x=kx, y=ky, z=kz, mode="lines", hoverinfo="skip")
p2 = scatter3d(; x=kx[adc_onmask], y=ky[adc_onmask], z=kz[adc_onmask], mode="markers")
p3 = scatter3d(; x=[0.], y=[0.], z=[0.], marker=attr(symbol="cross",size=10,color="red"))
display(plot([p1; p2; p3]))

# Simulate for new simple simulation
t, Δt, rfa, rfΔf, gxa, gya, gza, rf_onmask, gx_onmask, gy_onmask, gz_onmask, adc_onmask, adct, magxy, magz, sig = komasim(seq, obj)
display(plot([scatter(;x=t, y=abs.([sum(magxy; dims=1)...]), mode="lines+markers", name="sum(magxy)"); scatter(;x=t, y=abs.(sig), mode="lines+markers", name="sig")]))

# Simulate for old simulation
simParams = KomaMRICore.default_sim_params()
simParams["return_type"] = "mat"
told = KomaMRICore.get_adc_sampling_times(seq)
sigold = simulate(obj, seq, sys; simParams)
plot([scatter(;x=told, y=abs.(sigold[:,1,1]), mode="lines+markers", name="")])

# Compare new and old simulations
display(plot([scatter(;x=t, y=abs.(sig), mode="lines+markers"); scatter(;x=told, y=abs.(sigold[:,1,1]), mode="lines+markers")]))

using KomaMRI
using MAT

seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq")
seq = read_seq(seq_file)
sys = Scanner()
obj = brain_phantom2D()

############################################################################################
### Start Scanner

sys_dict = Dict("B0" => sys.B0,
                "B1" => sys.B1,
                "Gmax" => sys.Gmax,
                "Smax" => sys.Smax,
                "ADC_dt" => sys.ADC_Δt,
                "seq_dt" => sys.seq_Δt,
                "GR_dt" => sys.GR_Δt,
                "RF_dt" => sys.RF_Δt,
                "RF_ring_down_T" => sys.RF_ring_down_T,
                "RF_dead_time_T" => sys.RF_dead_time_T,
                "ADC_dead_time_T" => sys.ADC_dead_time_T)

### End Scanner
############################################################################################

############################################################################################
### Start Phantom

obj_dict = Dict("name" => obj.name,
                "columns" => ["x", "y", "z", "rho", "T1", "T2", "T2s", "delta_omega"],
                "data" => hcat(obj.x, obj.y, obj.z, obj.ρ, obj.T1, obj.T2, obj.T2s, obj.Δw))

### End Phantom
############################################################################################

############################################################################################
### Start Sequence

max_rf_samples=100
N = length(seq)
#O = size(seq.RF,1)
ΔT = KomaMRICore.durs(seq)
T0 = cumsum([0; ΔT],dims=1)
off_val = Inf #This removes the unnecessary points in the plot

#GRADS
t1x = vcat([KomaMRICore.get_theo_t(seq.GR[1,i]) .+ T0[i] for i=1:N]...)
t1y = vcat([KomaMRICore.get_theo_t(seq.GR[2,i]) .+ T0[i] for i=1:N]...)
t1z = vcat([KomaMRICore.get_theo_t(seq.GR[3,i]) .+ T0[i] for i=1:N]...)
Gx =  vcat([KomaMRICore.get_theo_A(seq.GR[1,i];off_val) for i=1:N]...)
Gy =  vcat([KomaMRICore.get_theo_A(seq.GR[2,i];off_val) for i=1:N]...)
Gz =  vcat([KomaMRICore.get_theo_A(seq.GR[3,i];off_val) for i=1:N]...)
GRADS = hcat(t1x, t1y, t1z, Gx, Gy, Gz)
#RFS
t2 =  vcat([KomaMRICore.get_theo_t(seq.RF[1,i];max_rf_samples) .+ T0[i] for i=1:N]...)
R =   vcat([KomaMRICore.get_theo_A(r;off_val,max_rf_samples) for r = seq.RF]...)
RFS = hcat(t2, R)
#ADC
t3 =  vcat([KomaMRICore.get_theo_t(seq.ADC[i])  .+ T0[i] for i=1:N]...)
D =   vcat([KomaMRICore.get_theo_A(d;off_val) for d = seq.ADC]...)
ADCS = hcat(t3, D)

seq_dict = Dict("GRAD" => GRADS,
                "RF" => RFS,
                "ADC" => ADCS,
                "DUR" => seq.DUR,
                "DEF" => seq.DEF)

### End Sequence
############################################################################################

############################################################################################
### Start kspace

kspace, kspace_adc = get_kspace(seq; Δt=1)

### End kspace
############################################################################################

############################################################################################
### Start M0

#Times
dt = 1
t, Δt = KomaMRICore.get_uniform_times(seq, dt)
#kx,ky
ts = t .+ Δt
k, _ =  KomaMRICore.get_kspace(seq; Δt=dt)

momentum0 = hcat(t, k)

### End M0
############################################################################################

############################################################################################
### Start simulation

raw = simulate(obj, seq, sys)

simParams = raw.params["userParameters"]

not_Koma = raw.params["systemVendor"] != "KomaMRI.jl"
t = Float64[]
signal = ComplexF64[]
current_t0 = 0
for p in raw.profiles
	dt = p.head.sample_time_us != 0 ? p.head.sample_time_us * 1e-3 : 1
	t0 = p.head.acquisition_time_stamp * 1e-3 #This parameter is used in Koma to store the time offset
	N =  p.head.number_of_samples != 0 ? p.head.number_of_samples : 1
	if not_Koma
		t0 = current_t0 * dt
		current_t0 += N
	end
	if N != 1
		append!(t, t0.+(0:dt:dt*(N-1)))
	else
		append!(t, t0)
	end
	append!(signal, p.data[:,1]) #Just one coil
	#To generate gap
	append!(t, t[end])
	append!(signal, [Inf + Inf*1im])
end

raw_dict = hcat(t, signal)

### End Simulation
############################################################################################

############################################################################################
### Start Reconstruction

# Get the acquisition data
acq = AcquisitionData(raw)
acq.traj[1].circular = false
Nx, Ny = raw.params["reconSize"][1:2]
reconParams = Dict{Symbol,Any}(:reco=>"direct", :reconSize=>(Nx, Ny))
rec = reconstruction(acq, reconParams)
image = rec.data[:,:,:]

reconParams_dict = Dict("reco" => reconParams[:reco],
                        "Nx" => reconParams[:reconSize][1],
                        "Ny" => reconParams[:reconSize][2])

### End Reconstruction
############################################################################################

matwrite("test.mat", Dict("scanner" => sys_dict,
                          "phantom" => obj_dict,
                          "sequence" => seq_dict,
                          "kspace" => kspace,
                          "kspace_adc" => kspace_adc,
                          "momentum0" => momentum0,
                          "sim_params" => simParams,
                          "raw" => raw_dict,
                          "rec_params" => reconParams_dict,
                          "image" => image))

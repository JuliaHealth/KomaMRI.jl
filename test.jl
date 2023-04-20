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

################################################
# Simple test
B1=1
N=6
T=1e-6
seq = Sequence()
seq += ADC(N, T)
seq += ADC(N, T)
for i in 1:2
    seq += RF(B1,T)
    seq += ADC(N, T)
    #seq += Delay(100*T)
end
seq += Grad(B1,T)

sys = Scanner()
obj = Phantom{Float64}(x=[0],T1=[T],T2=[T])


raw = simulate(obj, seq, sys)

plot_signal(raw)




####
using Blink
using Interact

global matfolder = pwd()
global mat_obs = Observable{Any}(matfolder)

mergeclasses(args...) = join(args, ' ')

customopendialog(; value = String[], label = "Open", icon = "far fa-folder-open", kwargs...) =
    customdialog(js"showOpenDialog"; value = value, label = label, icon = icon, kwargs...)

function customdialog(dialogtype; value, theme = gettheme(), className = "", label = "dialog", icon = nothing, options...)
    (value isa AbstractObservable) || (value = Observable(value))
    scp = Scope()
    setobservable!(scp, "output", value)
    clicks = Observable(scp, "clicks", 0)
    callback = @js function (val)
        $value[] = val
    end
    onimport(scp, js"""
    function () {
        const { dialog } = require('electron').remote;
        this.dialog = dialog;
    }
    """)
    onjs(clicks, js"""
    function (val) {
        console.log(this.dialog.$dialogtype($options, $callback));
    }
    """)
    className = mergeclasses("is-medium button", className)
    content = if icon === nothing
        (label,)
    else
        iconNode = node(:span, node(:i, className = icon), className = "icon")
        (iconNode, node(:span, label))
    end
    btn = node(:button, content...,
        events=Dict("click" => @js event -> ($clicks[] = $clicks[] + 1)),
        className = className)
    scp.dom = btn
    slap_design!(scp, theme)
    Widget{:dialog}([]; output = value, scope = scp, layout = Widgets.scope)
end

loadbutton = filepicker(; properties = ["openDirectory"])
dialogbutton = customopendialog(; properties = ["openDirectory"])
ui = vbox( # put things one on top of the other
    loadbutton,
    dialogbutton
)

#map!(f->if f!="" #Assigning function of data when load button (opendialog) is changed
#            println(f)
#            global matfolder = f[1]
#            @js_ w (@var name = $(basename(f[1]));
#            document.getElementById("folname").innerHTML=name)
#            matfolder
#        else
#            matfolder #default sequence
#        end
#    , mat_obs, dialogbutton)

w = Window()
body!(w, ui);



###############
seq = Sequence()  # empty sequence
seq += exc        # adding RF-only block
seq += acq        # adding ADC-only block
seq += Sequence([Grad(0, 0.1), Grad(0, 0.1), Grad(1, 0.1)])
p1 = plot_seq(seq; slider=false, height=300)

####################
using KomaMRI, Blink, Interact

darkmode = true
img_obs = Observable{Any}(image)
image = [0.0im; 2;; 3; 4;; 5; 6;;; 7; 8;; 9; 10;; 11; 12]
#image = [0.0im; 0.0im;; 0.0im; 0.0im;; 0.0im; 0.0im;;; 0.0im; 0.0im;; 0.0im; 0.0im;; 0.0im; 0.0im]

function plotInteract()
    @manipulate for slice = 1:size(image,3)
        aux = abs.(image) * prod(size(image)[1:2])
        plot_image(aux[:,:,slice],zmin=minimum(aux[:]),zmax=maximum(aux[:]);darkmode,title="Reconstruction ($slice/$(size(image,3)))")
    end
end

plt = Observable{Any}(plotInteract())

map!(t-> plotInteract(), plt, img_obs)
ui = dom"div"(plt)

w = Window()
body!(w, ui)
##################
using KomaMRI, Blink, Interact, PlotlyJS

sys = Scanner()
B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
N = 101
FOV = 23e-2
EPI = PulseDesigner.EPI(FOV, N, sys)
TE = 30e-3
d1 = TE-dur(EPI)/2-dur(EX)
if d1 > 0 DELAY = Delay(d1) end
seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI

darkmode = true
phantom = brain_phantom3D()
pha_obs = Observable{Phantom}(phantom)

loadbutton = filepicker()
function aux(data)
    brain_phantom2D()
end

map!(aux, pha_obs, loadbutton)


columnbuttons = Observable{Any}(dom"div"())
plt = Observable{Any}(PlotlyJS.plot());

function plotInteract(ph, key)
    @manipulate for t0_ms = range(0,dur(seq),5)*1e3
        t0_ms = 0.0
        plot_phantom_map(ph, key; t0=t0_ms, darkmode)
    end
end

function makebuttons(ph)
    prop = propertynames(ph)[5:end-3] #Removes name,x,y,z and ux,uy,uz
    propnm = [s for s=string.(prop)]
    buttons = button.(propnm)
    for (btn, key) in zip(reverse(buttons), reverse(prop))
        map!(t -> begin plotInteract(ph, key) end, plt, btn)
    end
    dom"div"(hbox(buttons))
end



map!(makebuttons, columnbuttons, pha_obs)

w = Window()
ui = dom"div"(vbox(loadbutton, columnbuttons, plt));
body!(w, ui)
####################

using CSV, DataFrames, Interact, Plots
loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(DataFrame)
plt = Observable{Any}(plot())
function aux(data)
    CSV.read(data, DataFrame)
end

map!(aux, data, loadbutton)

function makebuttons(df)
    buttons = button.(string.(names(df)))
    for (btn, name) in zip(buttons, names(df))
        map!(t -> histogram(df[!, name]), plt, btn)
    end
    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)

ui = dom"div"(loadbutton, columnbuttons, plt)

w = Window()
body!(w, ui)

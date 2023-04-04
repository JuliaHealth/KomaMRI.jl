function KomaUI(;dark=true,frame=true, phantom_mode="2D", sim=Dict{String,Any}(), rec=Dict{Symbol,Any}(), devTools=false)
## ASSETS
path = @__DIR__
assets = AssetRegistry.register(dirname(path*"/ui/assets/"))
scripts = AssetRegistry.register(dirname(path*"/ui/scripts/"))
css = AssetRegistry.register(dirname(path*"/ui/css/"))
# Assets
background = assets*"/spiral-bg.svg" #In Windows joinpath causes problems "/assetserver/...-assets\Logo.png"
logo = joinpath(assets, "Logo_dark.svg")
icon = joinpath(assets, "Icon.svg")
loading = joinpath(assets, "Loading.gif")
# JS
bsjs = joinpath(scripts, "bootstrap.bundle.min.js") #this already has Popper
bscss = joinpath(css,"bootstrap.min.css")
bsiconcss = joinpath(css,"bootstrap-icons.css")
jquery = joinpath(scripts,"jquery-3.4.1.slim.min.js")
# mathjaxsetup = joinpath(scripts, "mathjaxsetup.js")
# KaTeX
katexrender = joinpath(scripts, "auto-render.min.js")
katexjs = joinpath(scripts,"katex.min.js")
katexcss = joinpath(css,"katex.min.css")
# User defined JS and CSS
customcss = joinpath(css,"custom.css")
customjs = joinpath(scripts,"custom.js")
customjs2 = joinpath(scripts,"custom2.js")
sidebarcss = joinpath(css,"sidebars.css")
# Custom icons
icons = joinpath(css,"icons.css")
## WINDOW
global w = Blink.Window(Dict(
    "title"=>"KomaUI",
    "autoHideMenuBar"=>true,
    "frame"=>frame, #removes title bar
    "node-integration" => true,
    :icon=>path*"/ui/assets/Logo_icon.png",
    "width"=>1200,
    "height"=>800,
    "webPreferences" => Dict("devTools" => devTools)
    ),async=false);
## LOADING BAR
buffericon = """<div class="spinner-border spinner-border-sm text-light" role="status"></div>"""
progressbar = """
<div class="progress" style="background-color: #27292d;">
  <div id="simul_progress" class="progress-bar" role="progressbar" style="width: 0%; transition:none;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">0%</div>
</div>
"""
## NAV BAR
sidebar = open(f->read(f, String), path*"/ui/html/sidebar.html")
sidebar = replace(sidebar, "LOGO"=>logo)
## CONTENT
index = open(f->read(f, String), path*"/ui/html/index.html")
index = replace(index, "ICON"=>icon)
#index = replace(index, "BACKGROUND_IMAGE"=>background)
## CSS
loadcss!(w, bscss)
loadcss!(w, bsiconcss)
loadcss!(w, customcss)
loadcss!(w, sidebarcss)
# KATEX
loadcss!(w, katexcss)
loadjs!(w, katexjs)
loadjs!(w, katexrender)
# JQUERY, BOOSTRAP JS
loadjs!(w, customjs)    #must be before jquery
loadjs!(w, jquery)
loadjs!(w, bsjs)        #after jquery
loadjs!(w, customjs2)   #must be after jquery
# LOAD ICONS
loadcss!(w, icons)
## INITAL VARIABLES
#PHANTOM init
@info "Loading Phantom (default)"
if phantom_mode == "3D"
    global phantom = brain_phantom3D()
else
    global phantom = brain_phantom2D()
end
phantom.Δw .*= 0
println("Phantom object \"$(phantom.name)\" successfully loaded!")
#SCANNER init
@info "Loading Scanner (default)"
global sys = Scanner()
println("B0 = $(sys.B0) T")
println("Gmax = $(round(sys.Gmax*1e3,digits=2)) mT/m")
println("Smax = $(sys.Smax) mT/m/ms")
#SEQ init
@info "Loading Sequence (default) "
B1 = sys.B1; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
N = 101
FOV = 23e-2
EPI = PulseDesigner.EPI(FOV, N, sys)
TE = 30e-3
d1 = TE-dur(EPI)/2-dur(EX)
if d1 > 0 DELAY = Delay(d1) end
global seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI
seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
#Init
global darkmode = dark
global raw_ismrmrd = RawAcquisitionData(Dict(
    "systemVendor" => "",
    "encodedSize" => [2,2,1],
    "reconSize" => [2,2,1],
    "number_of_samples" => 4,
    "encodedFOV" => [100.,100.,1],
    "trajectory" => "other"),
    [KomaMRICore.Profile(AcquisitionHeader(trajectory_dimensions=2, sample_time_us=1),
        [0. 0. 1 1; 0 1 1 1]./2, [0.; 0im; 0; 0;;])])
global rawfile = ""
global image =  [0.0im 0.; 0. 0.]
global kspace = [0.0im 0.; 0. 0.]
#Reco
default = Dict{Symbol,Any}(:reco=>"direct") #, :iterations=>10, :λ=>1e-5,:solver=>"admm",:regularization=>"TV")
global recParams = merge(default, rec)
#Simulation
default = Dict{String,Any}()
global simParams = merge(default, sim)
#GPUs
@info "Loading GPUs"
KomaMRICore.print_gpus()
#OBERSVABLES
global seq_obs = Observable{Sequence}(seq)
global pha_obs = Observable{Phantom}(phantom)
global sig_obs = Observable{RawAcquisitionData}(raw_ismrmrd)
global img_obs = Observable{Any}(image)

function export2matsequence()
	max_rf_samples=100
    N = length(seq)
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
    matwrite("sequence.mat", Dict("sequence" => seq_dict))
end

function export2matkspace()
    kspace, kspace_adc = get_kspace(seq; Δt=1)
    matwrite("kspace.mat", Dict("kspace" => kspace, "kspace_adc" => kspace_adc))
end

function export2matmoment0()
    dt = 1
    t, Δt = KomaMRICore.get_uniform_times(seq, dt)
    ts = t .+ Δt
    k, _ =  KomaMRICore.get_kspace(seq; Δt=dt)
    moment0 = hcat(t, k)
    matwrite("moment0.mat", Dict("moment0" => moment0))
end

function export2matphantom()
    phantom_dict = Dict("name" => phantom.name,
                "columns" => ["x", "y", "z", "rho", "T1", "T2", "T2s", "delta_omega"],
                "data" => hcat(phantom.x, phantom.y, phantom.z, phantom.ρ, phantom.T1, phantom.T2, phantom.T2s, phantom.Δw))
    matwrite("phantom.mat", Dict("phantom" => phantom_dict))
end

function export2matscanner()
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
    matwrite("scanner.mat", Dict("scanner" => sys_dict))
end

function export2matraw()
    if haskey(raw_ismrmrd.params, "userParameters")
        matwrite("sim_params.mat", Dict("sim_params" => raw_ismrmrd.params["userParameters"]))

        not_Koma = raw_ismrmrd.params["systemVendor"] != "KomaMRI.jl"
        t = Float64[]
        signal = ComplexF64[]
        current_t0 = 0
        for p in raw_ismrmrd.profiles
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
        matwrite("raw.mat", Dict("raw" => raw_dict))
    end

end

function export2matimage()
    if haskey(recParams, :reconSize)
        recParams_dict = Dict("reco" => recParams[:reco],
                            "Nx" => recParams[:reconSize][1],
                            "Ny" => recParams[:reconSize][2])
        matwrite("rec_params.mat", Dict("rec_params" => recParams_dict))
    end

    matwrite("image.mat", Dict("image" => image))
end

function export2mat(;type="all")
    if type=="all"
        export2matsequence()
        export2matkspace()
        export2matmoment0()
        export2matphantom()
        export2matscanner()
    elseif type=="sequence"
		export2matsequence()
        export2matkspace()
        export2matmoment0()
	elseif type=="phantom"
		export2matphantom()
    elseif type=="scanner"
		export2matscanner()
    elseif type=="raw"
		export2matraw()
    elseif type=="image"
		export2matimage()
	end
end

handle(w, "matall") do args...
    export2mat(;type="all")
end
handle(w, "matsequence") do args...
    export2mat(;type="sequence")
end
handle(w, "matphantom") do args...
    export2mat(;type="phantom")
end
handle(w, "matscanner") do args...
    export2mat(;type="scanner")
end
handle(w, "matraw") do args...
    export2mat(;type="raw")
end
handle(w, "matimage") do args...
    export2mat(;type="image")
end

## MENU FUNCTIONS
handle(w, "index") do args...
    content!(w, "div#content", index)
end
handle(w, "pulses_seq") do args...
    include(path*"/ui/PulsesGUI_seq.jl")
end
handle(w, "pulses_kspace") do args...
    include(path*"/ui/PulsesGUI_kspace.jl")
end
handle(w, "pulses_M0") do args...
    include(path*"/ui/PulsesGUI_M0.jl")
end
handle(w, "phantom") do args...
    include(path*"/ui/PhantomGUI.jl")
end
handle(w, "sig") do args...
    include(path*"/ui/SignalGUI.jl")
end
handle(w, "reconstruction_absI") do args...
    include(path*"/ui/ReconGUI_absI.jl")
end
handle(w, "reconstruction_angI") do args...
    include(path*"/ui/ReconGUI_angI.jl")
end
handle(w, "reconstruction_absK") do args...
    include(path*"/ui/ReconGUI_absK.jl")
end
handle(w, "sim_params") do args...
    include(path*"/ui/SimParams_view.jl")
end
handle(w, "rec_params") do args...
    include(path*"/ui/RecParams_view.jl")
end
handle(w, "simulate") do args...
    # @js_ w document.getElementById("simulate!").prop("disabled", true); #Disable button during SIMULATION
    @js_ w (@var progressbar = $progressbar; document.getElementById("simulate!").innerHTML=progressbar)
    #To SequenceGUI
    global raw_ismrmrd = simulate(phantom, seq, sys; simParams, w)
    #After simulation go to RECON
    @js_ w document.getElementById("simulate!").innerHTML="Simulate!"
    #EXPORT to ISMRMRD -> To SignalGUI
    global rawfile = tempdir()*"/Koma_signal.mrd"
    @info "Exporting to ISMRMRD file: $rawfile"
    global sig_obs[] = raw_ismrmrd
    fout = ISMRMRDFile(rawfile)
    KomaMRICore.save(fout, raw_ismrmrd)
    #Message
    sim_time = raw_ismrmrd.params["userParameters"]["sim_time_sec"]
    @js_ w (@var sim_time = $sim_time;
    @var name = $(phantom.name);
    document.getElementById("rawname").innerHTML=name;
    Toasty("1", """Simulation successfull<br>Time: <a id="sim_time"></a> s""" ,"""
    <ul>
        <li>
            <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('sig', 1)"><i class="fa fa-search"></i></button>
            Updating <b>Raw signal</b> plots ...
        </li>
        <li>
            <button class="btn btn-primary btn-circle btn-circle-sm m-1" onclick="Blink.msg('recon', 1)"><i class="bi bi-caret-right-fill"></i></button>
            Ready to <b>reconstruct</b>?
        </li>
    </ul>
    """);
    document.getElementById("sim_time").innerHTML=sim_time;
    )
    # @js_ w document.getElementById("simulate!").prop("disabled", false); #Re-enable button
    # @js_ w (@var button = document.getElementById("recon!"); @var bsButton = @new bootstrap.Button(button); vsButton.toggle())
end
handle(w, "recon") do args...
    # Update loading icon for button
    @js_ w (@var buffericon = $buffericon; document.getElementById("recon!").innerHTML=buffericon)
    #IMPORT ISMRMRD raw data
    raw_ismrmrd.profiles = raw_ismrmrd.profiles[getproperty.(getproperty.(raw_ismrmrd.profiles, :head), :flags) .!= 268435456] #Extra profile in JEMRIS simulations
    acqData = AcquisitionData(raw_ismrmrd)
    acqData.traj[1].circular = false #Removing circular window
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ maximum(2*abs.(acqData.traj[1].nodes[:])) #Normalize k-space to -.5 to .5 for NUFFT
    Nx, Ny = raw_ismrmrd.params["reconSize"][1:2]
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    #Reconstruction
    @info "Running reconstruction of $rawfile"
    aux = @timed reconstruction(acqData, recParams)
    global image  = reshape(aux.value.data,Nx,Ny,:)
    global kspace = fftc(reshape(aux.value.data,Nx,Ny,:))
    # global img_obs[] = image
    #After Recon go to Image
    recon_time = aux.time
    @js_ w document.getElementById("recon!").innerHTML="Reconstruct!"
    @js_ w (@var recon_time = $recon_time;
    @var name = $(phantom.name);
    document.getElementById("imaname").innerHTML=name;
    Toasty("2", """Reconstruction successfull<br>Time: <a id="recon_time"></a> s""" ,"""
    <ul>
        <li>
            <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('reconstruction_absI', 1)"><i class="fa fa-search"></i></button>
            Updating <b>Reconstruction</b> plots ...
        </li>
    </ul>
    """
    );
    document.getElementById("recon_time").innerHTML=recon_time;
    )
end
handle(w, "close") do args...
    global darkmode = nothing

    global phantom = nothing
    global seq = nothing
    global sys = nothing

    global raw_ismrmrd = nothing
    global rawfile = nothing

    global image = nothing
    global kspace = nothing

    global simParams = nothing
    global recParams = nothing

    global seq_obs = nothing
    global pha_obs = nothing
    global sig_obs = nothing
    global img_obs = nothing
    close(w)
end
#Update GUI's home
w = body!(w, *(sidebar,index), async=false)
if darkmode
    @js_ w document.getElementById("main").style="background-color:rgb(13,16,17);"
end
#Sequence observable
load_seq = filepicker(".seq (Pulseq)/.seqk (Koma)"; accept=".seq,.seqk")
map!(f->if f!="" #Assigning function of data when load button (filepicker) is changed
            if splitext(f)[end]==".seqk" #Koma
                global seq = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"seq")
            elseif splitext(f)[end]==".seq" #Pulseq
                global seq = read_seq(f) #Pulseq read
            end
            @js_ w (@var name = $(basename(f));
            document.getElementById("seqname").innerHTML=name;
            Toasty("0", "Loaded <b>"+name+"</b> successfully", """
            <ul>
                <li>
                    <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('pulses_seq', 1)"><i class="fa fa-search"></i></button>
                    Updating <b>Sequence</b> plots ...
                </li>
                <li>
                    <button class="btn btn-primary btn-circle btn-circle-sm m-1" onclick="Blink.msg('simulate', 1)"><i class="bi bi-caret-right-fill"></i></button>
                    Ready to <b>simulate</b>?
                </li>
            </ul>
            """))
            seq
        else
            seq #default sequence
        end
    , seq_obs, load_seq)
w = content!(w, "#seqfilepicker", load_seq, async=false)
#Phantom observable
load_pha = filepicker(".phantom (Koma)/.h5 (JEMRIS)"; accept=".phantom,.h5")
map!(f->if f!="" #Assigning function of data when load button (filepicker) is changed
            if splitext(f)[end]==".phantom" #Koma
                global phantom = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"phantom")
            elseif splitext(f)[end]==".h5" #JEMRIS
                global phantom = read_phantom_jemris(f)
            end
            @js_ w (@var name = $(basename(f));
            document.getElementById("phaname").innerHTML=name;
            Toasty("0", "Loaded <b>"+name+"</b> successfully", """
            <ul>
                <li>
                    <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('phantom', 1)"><i class="fa fa-search"></i></button>
                    Updating <b>Phantom</b> plots ...
                </li>
                <li>
                    <button class="btn btn-primary btn-circle btn-circle-sm m-1" onclick="Blink.msg('simulate', 1)"><i class="bi bi-caret-right-fill"></i></button>
                    Ready to <b>simulate</b>?
                </li>
            </ul>
            """))
            phantom
        else
            phantom #default
        end
    , pha_obs, load_pha)
w = content!(w, "#phafilepicker", load_pha, async=false)
#Signal observable
load_sig = filepicker(".h5/.mrd (ISMRMRD)"; accept=".h5,.mrd")
map!(f->if f!="" #Assigning function of data when load button (filepicker) is changed
            fraw = ISMRMRDFile(f)
            global rawfile = f
            global raw_ismrmrd = RawAcquisitionData(fraw)

            not_Koma = raw_ismrmrd.params["systemVendor"] != "KomaMRI.jl"
            if not_Koma
                @warn "ISMRMRD files generated externally could cause problems during the reconstruction. We are currently improving compatibility."
            end

            @js_ w (@var name = $(basename(f));
            document.getElementById("rawname").innerHTML=name;
            Toasty("0", "Loaded <b>"+name+"</b> successfully", """
            <ul>
                <li>
                    <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('sig', 1)"><i class="fa fa-search"></i></button>
                    Updating <b>Raw data</b> plots ...
                </li>
                <li>
                    <button class="btn btn-success btn-circle btn-circle-sm m-1" onclick="Blink.msg('recon', 1)"><i class="bi bi-caret-right-fill"></i></button>
                    Ready to <b>reconstruct</b>?
                </li>
            </ul>
            """))
            raw_ismrmrd
        else
            raw_ismrmrd #default
        end
    , sig_obs, load_sig)
w = content!(w, "#sigfilepicker", load_sig, async=false)
#Update Koma version
version = string(KomaMRICore.__VERSION__)
content!(w, "#version", version, async=false)
@info "Currently using KomaMRICore v$version"

content!(w, "#phaname", phantom.name, async=false)

return w
end

"""Updates KomaUI's simulation progress bar."""
function update_blink_window_progress!(w::Window, block, Nblocks)
    progress = string(floor(Int, block / Nblocks * 100))
    @js_ w (@var progress = $progress;
    document.getElementById("simul_progress").style.width = progress + "%";
    document.getElementById("simul_progress").innerHTML = progress + "%";
    document.getElementById("simul_progress").setAttribute("aria-valuenow", progress))
    return nothing
end

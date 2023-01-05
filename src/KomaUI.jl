function KomaUI(;dark=true,frame=true, phantom_mode="2D", sim=Dict{String,Any}(), rec=Dict{Symbol,Any}())
## ASSETS
path = @__DIR__
assets = AssetRegistry.register(dirname(path*"/ui/assets/"))
scripts = AssetRegistry.register(dirname(path*"/ui/scripts/"))
css = AssetRegistry.register(dirname(path*"/ui/css/"))
# Assets
background = assets*"/spiral-bg.svg" #In Windows joinpath causes problems "/assetserver/...-assets\Logo.png"
logo = joinpath(assets, "Logo_dark.svg")
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
    "height"=>800
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
index = replace(index, "BACKGROUND_IMAGE"=>background)
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
    [KomaMRI.Profile(AcquisitionHeader(trajectory_dimensions=2, sample_time_us=1),
        [0. 0. 1 1; 0 1 1 1]./2, [0.; 0im; 0; 0;;])])
global rawfile = ""
global image =  [0.0im 0.; 0. 0.]
global kspace = [0.0im 0.; 0. 0.]
#Reco
default = Dict{Symbol,Any}(:reco=>"direct") #, :iterations=>10, :λ=>1e-5,:solver=>"admm",:regularization=>"TV")
global recParams = merge(default, rec)
#Simulation
default = Dict{String,Any}("gpu"=>has_cuda(), "gpu_device"=>0, "Nthreads"=>has_cuda() ? 1 : Threads.nthreads())
global simParams = merge(default, sim)
#GPUs
if simParams["gpu"]
    @info "Loading GPUs"
    print_gpus()
    device!(simParams["gpu_device"]) #By default it uses first GPU, multiGPU not supported yet
end
#OBERSVABLES
global seq_obs = Observable{Sequence}(seq)
global pha_obs = Observable{Phantom}(phantom)
global sig_obs = Observable{RawAcquisitionData}(raw_ismrmrd)
global img_obs = Observable{Any}(image)
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
    save(fout, raw_ismrmrd)
    #Message
    sim_time = raw_ismrmrd.params["userParameters"]["sim_time_sec"]
    @js_ w (@var sim_time = $sim_time; 
    Toasty("1", """Simulation successfull<br>Time: <a id="sim_time"></a> s""" ,"""
    <ul>
        <li>
            <button class="btn btn-dark btn-circle btn-circle-sm m-1" onclick="Blink.msg('sig', 1)"><i class="fa fa-search"></i></button>
            Updating <b>Raw signal</b> plots ...
        </li>
        <li>
            <button class="btn btn-success btn-circle btn-circle-sm m-1" onclick="Blink.msg('recon', 1)"><i class="bi bi-caret-right-fill"></i></button>
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
    aux = @timed MRIReco.reconstruction(acqData, recParams)
    global image  = reshape(aux.value.data,Nx,Ny,:)
    global kspace = fftc(reshape(aux.value.data,Nx,Ny,:))
    # global img_obs[] = image
    #After Recon go to Image
    recon_time = aux.time
    @js_ w document.getElementById("recon!").innerHTML="Reconstruct!"
    @js_ w (@var recon_time = $recon_time; 
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
version =  string(VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"]))
content!(w, "#version", version, async=false)
@info "Currently using KomaMRI v$version"

nothing
end

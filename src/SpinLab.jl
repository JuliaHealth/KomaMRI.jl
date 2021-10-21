function SpinLab(;frame=true)
## ASSETS
path = @__DIR__
assets = AssetRegistry.register(dirname(path*"/ui/assets/"))
scripts = AssetRegistry.register(dirname(path*"/ui/scripts/"))
css = AssetRegistry.register(dirname(path*"/ui/css/"))

logo = joinpath(assets, "Logo.png")
loading = joinpath(assets, "Loading.gif")
popper = joinpath(scripts, "popper.min.js")
bsjs = joinpath(scripts, "bootstrap.min.js")
bscss = joinpath(css,"bootstrap.min.css")
jquery = joinpath(scripts,"jquery-3.4.1.slim.min.js")
# KaTeX
katexrender = joinpath(scripts, "auto-render.min.js")
katexjs = joinpath(scripts,"katex.min.js")
katexcss = joinpath(css,"katex.min.css")
# User defined
customcss = joinpath(css,"custom.css")
customjs = joinpath(scripts,"custom.js")
customjs2 = joinpath(scripts,"custom2.js")
# Custom icons
icons = joinpath(css,"icons.css")
# Others
imphantom = joinpath(assets, "phantom.png")
imscanner = joinpath(assets, "scanner.png")
impulses = joinpath(assets, "pulses.png")
## WINDOW
global w = Blink.Window(Dict(
    "title"=>"SpinLab",
    # "autoHideMenuBar"=>false,
    "frame"=>frame, #removes title bar
    "node-integration" => true,
    :icon=>path*"/ui/assets/Logo_icon.png"
    ),async=false);
## LOADING BAR
loadbar = """<center><img src="$loading" width="100px" class="align-middle"></center><br>"""
## NAV BAR
navbar = open(f->read(f, String), path*"/ui/html/navbar.html")
navbar = replace(navbar, "LOGO"=>logo)
## CONTENT
index = open(f->read(f, String), path*"/ui/html/index.html")
## FOOTER
footer = open(f->read(f, String), path*"/ui/html/footer.html")
## CSS
loadcss!(w, bscss)
loadcss!(w, customcss)
# KATEX
loadcss!(w, katexcss)
loadjs!(w, katexjs)
loadjs!(w, katexrender)
# JQUERY, BOOSTRAP JS
# mathjax = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-AMS-MML_SVG"
# loadjs!(w, mathjax)
loadjs!(w, customjs)    #must be before jquery
loadjs!(w, jquery)
loadjs!(w, customjs2)   #must be after jquery
loadjs!(w, popper)
loadjs!(w, bsjs)        #after jquery
# LOAD IMAGES AND ICONS
loadcss!(w, icons)
index = replace(index, "PHANTOM"=>imphantom)
index = replace(index, "SCANNER"=>imscanner)
index = replace(index, "PULSES"=>impulses)
## MENU FUNCTIONS
handle(w, "index") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    body!(w,*(navbar,index,footer))
end
handle(w, "pulses") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include(path*"/ui/PulsesGUI.jl")
end
handle(w, "phantom") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include(path*"/ui/PhantomGUI.jl")
end
handle(w, "sim") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include(path*"/ui/SimulatorGUI.jl")
end
handle(w, "recon") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include(path*"/ui/ReconGUI.jl")
end
handle(w, "simulate") do args...
    @info "Running simulation..."
    @js_ w (@var loading = $loadbar; document.getElementById("simulate").innerHTML=loading)
    Δt = 4e-6 #<- simulate param
    t = collect(0:Δt:MRIsim.dur(seq))
    Nphant, Nt = prod(size(phantom)), length(t)
    N_parts = floor(Int, Nphant*Nt/2.7e6)
    println("Dividing simulation in Nblocks=$N_parts")
    S = @time MRIsim.run_sim_time_iter(phantom,seq,t;N_parts)
    global signal = S ./prod(size(phantom)) #Acquired data
    #Recon, will be replaced by call to MRIReco.jl
    S = nothing #remove aux signal S
    Nx = Ny = 99 #hardcoded by now
    global kdata = reshape(signal,(Nx,Ny)) #Turning into kspace image
    global kdata[:,2:2:Ny] = kdata[Nx:-1:1,2:2:Ny] #Flip in freq-dir for the EPI
    global kdata = convert(Array{Complex{Float64},2},kdata)

    @js_ w document.getElementById("simulate").innerHTML="""<button type="button" onclick='Blink.msg("simulate", 1)' class="btn btn-primary btn-lg btn-block">Run simulation!</button>"""
    @js_ w Blink.msg("recon", 0)
end
handle(w, "close") do args...
    global phantom = nothing
    global seq = nothing
    global scanner = nothing

    close(w)
end
## PRINTING INFO
#PHANTOM init
@info "Loading Phantom (default)"
global phantom = brain_phantom2D(;axis="coronal")
println("Phantom object \"$(phantom.name)\" successfully loaded!")
#SEQ init
@info "Loading Sequence (default) "
B1 = 6e-6; durRF = π/2/(2π*γ*B1) #90-degree hard excitation pulse
EX = PulseDesigner.RF_hard(B1, durRF)
Gmax = 60e-3
EPI,_,_,_ = PulseDesigner.EPI_base(40/100, 99, 4e-6, Gmax)
TE = 25e-3 
d = delay(TE-dur(EPI)/2-dur(EX))
DELAY = Sequence([d;d])
global seq = EX + DELAY + EPI
println("EPI successfully loaded! (TE = $(TE*1e3) ms)")
#Init
global scanner = []
global signal = 0
global kdata = [0.0im 0.; 0. 0.]
#GPUs
if has_cuda()
    @info "Loading GPUs"
    print_gpus()
end
#Update GUI's home
body!(w,*(navbar,index,footer),async=false)
nothing
end
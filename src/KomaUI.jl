function KomaUI(;frame=true)
## ASSETS
path = @__DIR__
assets = AssetRegistry.register(dirname(path*"/ui/assets/"))
scripts = AssetRegistry.register(dirname(path*"/ui/scripts/"))
css = AssetRegistry.register(dirname(path*"/ui/css/"))

# Assets
imphantom = assets*"/phantom.png" #In Windows joinpath causes problems "/assetserver/...-assets\Logo.png"
imscanner = assets*"/scanner.png"
impulses = assets*"/pulses.png"
logo = joinpath(assets, "Logo_dark.png")
loading = joinpath(assets, "Loading.gif")
# JS 
popper = joinpath(scripts, "popper.min.js")
bsjs = joinpath(scripts, "bootstrap.min.js")
bscss = joinpath(css,"bootstrap.min.css")
jquery = joinpath(scripts,"jquery-3.4.1.slim.min.js")
# KaTeX
katexrender = joinpath(scripts, "auto-render.min.js")
katexjs = joinpath(scripts,"katex.min.js")
katexcss = joinpath(css,"katex.min.css")
# User defined JS and CSS
customcss = joinpath(css,"custom.css")
customjs = joinpath(scripts,"custom.js")
customjs2 = joinpath(scripts,"custom2.js")
# Custom icons
icons = joinpath(css,"icons.css")

## WINDOW
global w = Blink.Window(Dict(
    "title"=>"KomaUI",
    "autoHideMenuBar"=>true,
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
index = replace(index, "PHANTOM"=>imphantom)
index = replace(index, "SCANNER"=>imscanner)
index = replace(index, "PULSES"=>impulses)
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
# LOAD ICONS
loadcss!(w, icons)
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
handle(w, "sig") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include(path*"/ui/SignalGUI.jl")
end
handle(w, "recon") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include(path*"/ui/ReconGUI.jl")
end
handle(w, "simulate") do args...
    @info "Running simulation..."
    @js_ w (@var loading = $loadbar; document.getElementById("simulate").innerHTML=loading)
    #To SequenceGUI
    aux = simulate(phantom, seq, sys, simParams)
    #To SignalGUI
    global signal = aux[1]
    global t_interp = aux[2]
    #To ReconGUI
    global recParams = aux[3]
    @js_ w Blink.msg("recon", 0)
    @js_ w document.getElementById("simulate").innerHTML="""<button type="button" onclick='Blink.msg("simulate", 1)' class="btn btn-primary btn-lg btn-block">Run simulation!</button>"""
end
handle(w, "close") do args...
    global phantom = nothing
    global seq = nothing
    global scanner = nothing

    global signal = nothing
    global t_interp = nothing

    global simParams = nothing
    global recParams = nothing

    close(w)
end
## PRINTING INFO
#PHANTOM init
@info "Loading Phantom (default)"
global phantom = brain_phantom2D()
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
EPI = PulseDesigner.EPI(23e-2, 101, sys)
TE = 30e-3
d1 = TE-dur(EPI)/2-dur(EX)
if d1 > 0 DELAY = Delay(d1) end
global seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI
seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1,digits=4)*1e3
println("EPI successfully loaded! (TE = $(seq.DEF["TE"]) ms)")
#Init
global t_interp = 0
global signal = [0.0im 0. 0. 0.]
global image = [0.0im 0.; 0. 0.]
global simParams = Dict("step"=>"uniform","Δt"=>1e-3)
global recParams = Dict("Nx"=>2,"recon"=>"fft")
#GPUs
if has_cuda()
    @info "Loading GPUs"
    print_gpus()
end
#Update GUI's home
body!(w,*(navbar,index,footer),async=false)
nothing
end
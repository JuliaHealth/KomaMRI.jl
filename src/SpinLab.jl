function SpinLab(;frame=true)

## ASSETS
path = @__DIR__
logo = AssetRegistry.register(path*"/ui/assets/Logo.png")
loading = AssetRegistry.register(path*"/ui/assets/Loading.gif")
# jquery = AssetRegistry.register(path*"/ui/scripts/jquery-3.4.1.slim.min.js")
popper = AssetRegistry.register(path*"/ui/scripts/popper.min.js")
bsjs = AssetRegistry.register(path*"/ui/scripts/bootstrap.min.js")
bscss = AssetRegistry.register(path*"/ui/css/bootstrap.min.css")
katex = AssetRegistry.register(path*"/ui/scripts/auto-render.min.js")
customcss = AssetRegistry.register(path*"/ui/css/custom.css")
# customjs = AssetRegistry.register(path*"/ui/scripts/custom.js")
# customjstmp = AssetRegistry.register(path*"/ui/scripts/custom_tmp.js")

# custom icons
icons = AssetRegistry.register(path*"/ui/css/icons.css")
# iconstmp = AssetRegistry.register(path*"/ui/css/icons_tmp.css")
# fontseot = AssetRegistry.register(path*"/ui/css/fonts/icomoon.eot")
# fontsttf = AssetRegistry.register(path*"/ui/css/fonts/icomoon.ttf")
# fontswoff = AssetRegistry.register(path*"/ui/css/fonts/icomoon.woff")
# fontssvg = AssetRegistry.register(path*"/ui/css/fonts/icomoon.svg")
#others
imphantom = AssetRegistry.register(path*"/ui/assets/phantom.png")
imscanner = AssetRegistry.register(path*"/ui/assets/scanner.png")
impulses = AssetRegistry.register(path*"/ui/assets/pulses.png")
## WINDOW
global w = Blink.Window(Dict(
    "title"=>"SpinLab",
    # "darkTheme"=>true,
    # "autoHideMenuBar"=>false,
    "frame"=>frame, #removes title bar
    # "transparent"=>false,
    # "backgroundcolor"=>"#000000",
    # "node-integration" => true,
    # :show=>false,
    :icon=>path*"/ui/assets/Logo_icon.png",
    # "minHeight"=>700,
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
#KATEX
loadcss!(w,"https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.css")
loadjs!(w,"https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.js")
loadjs!(w, katex)
#JQUERY, BOOSTRAP JS
# customjs_tmp = open(f->read(f, String), customjs)
# customjs_tmp = replace(customjs_tmp, "JQUERY"=>path*"/ui/scripts/jquery-3.4.1.slim.min.js")
# open(path*"/ui/scripts/custom_tmp.js", "w") do f write(f, customjs) end
# loadjs!(w, customjs_tmp)
loadjs!(w, popper)
loadjs!(w, bsjs)
# loadjs_defer!(w, jquery)
#ICONS
# icons = open(f->read(f, String), path*"/ui/css/icons.css")
# icons = replace(icons, "fontseot"=>fontseot)
# icons = replace(icons, "fontsttf"=>fontsttf)
# icons = replace(icons, "fontswoff"=>fontswoff)
# icons = replace(icons, "fontssvg"=>fontssvg)

# navbar = replace(navbar, "fontseot"=>path*"/ui/css/fonts/icomoon.eot")
# navbar = replace(navbar, "fontsttf"=>path*"/ui/css/fonts/icomoon.ttf")
# navbar = replace(navbar, "fontswoff"=>path*"/ui/css/fonts/icomoon.woff")
# navbar = replace(navbar, "fontssvg"=>path*"/ui/css/fonts/icomoon.svg")
# open(path*"/ui/css/icons_tmp.css", "w") do f write(f, icons) end
loadcss!(w, "ui/css/icons.css")
# LOAD IMAGES
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
    t = collect(Δt:Δt:MRIsim.dur(seq))
    Nphant, Nt = prod(size(phantom)), length(t)
    N_parts = floor(Int, Nphant*Nt/2.7e6)
    println("Dividing simulation in Nblocks=$N_parts")
    S = @time MRIsim.run_sim2D_times_iter(phantom,seq,t;N_parts) #run_sim2D_times_iter run_sim2D_spin
    global signal = S[MRIsim.get_DAC_on(seq,t)]/prod(size(phantom)) #Acquired data
    #Recon
    S = nothing #remove aux signal S
    Nx = Ny = 99 #hardcoded by now
    global kdata = reshape(signal,(Nx,Ny)) #Turning into kspace image
    global kdata[:,2:2:Ny,:] = kdata[Nx:-1:1,2:2:Ny] #Flip in freq-dir for the EPI
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
Gmax = 60e-3
EPI,_,_,_ = PulseDesigner.EPI_base(40/100, 99, 4e-6, Gmax)
TE = 25e-3 
d = delay(TE-dur(EPI)/2)
DELAY = Sequence([d;d])
global seq = DELAY + EPI
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
#Update GUI
body!(w,*(navbar,index,footer),async=false)
nothing
end
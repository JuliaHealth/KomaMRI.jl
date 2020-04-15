function SpinLab(;frame=false)

## ASSETS
path = @__DIR__
logo = AssetRegistry.register(path*"/ui/assets/Logo.png")
loading = AssetRegistry.register(path*"/ui/assets/Loading.gif")
jquery = AssetRegistry.register(path*"/ui/scripts/jquery-3.4.1.slim.min.js")
popper = AssetRegistry.register(path*"/ui/scripts/popper.min.js")
bsjs = AssetRegistry.register(path*"/ui/scripts/bootstrap.min.js")
bscss = AssetRegistry.register(path*"/ui/css/bootstrap.min.css")
katex = AssetRegistry.register(path*"/ui/scripts/auto-render.min.js")
customcss = AssetRegistry.register(path*"/ui/css/custom.css")
customjs = AssetRegistry.register(path*"/ui/scripts/custom.js")
customjstmp = AssetRegistry.register(path*"/ui/scripts/custom_tmp.js")

# custom icons
icons = AssetRegistry.register(path*"/ui/css/icons.css")
iconstmp = AssetRegistry.register(path*"/ui/css/icons_tmp.css")
fontseot = AssetRegistry.register(path*"/ui/css/fonts/icomoon.eot")
fontsttf = AssetRegistry.register(path*"/ui/css/fonts/icomoon.ttf")
fontswoff = AssetRegistry.register(path*"/ui/css/fonts/icomoon.woff")
fontssvg = AssetRegistry.register(path*"/ui/css/fonts/icomoon.svg")
#others
imphantom = AssetRegistry.register(path*"/ui/assets/phantom.png")
## WINDOW
global w = Blink.Window(Dict(
    "title"=>"SpinLab",
    "darkTheme"=>true,
    "autoHideMenuBar"=>true,
    "frame"=>frame, #removes title bar
    "transparent"=>true,
    "backgroundColor"=>"#000",
    "node-integration" => true,
    :icon=>path*"/ui/assets/Logo_icon.png",
    "minHeight"=>500,
    # "minWidth"=>800
    ))
function loadjs_defer!(w, url)
  @js w @new Promise(function (resolve, reject)
    @var script = document.createElement("script")
    script.src = $url
    script.onload = resolve
    script.onerror = (e) -> reject(
                               Dict("name"=>"JSLoadError",
                                    "message"=>"failed to load " + this.src)
                                    )
    script.defer = true;
    document.head.appendChild(script)
  end)
end
## LOADING BAR
loadbar = """<center><img src="$loading" width="100px" class="align-middle"></center><br>"""
## NAV BAR
navbar = open(f->read(f, String), path*"/ui/html/navbar.html")
navbar = replace(navbar, "LOGO"=>logo)
## CONTENT
index = open(f->read(f, String), path*"/ui/html/index.html")
## FOOTER
footer = open(f->read(f, String), path*"/ui/html/footer.html")
## Update content
loadcss!(w, bscss)
loadcss!(w, customcss)
loadcss!(w,"https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.css")
loadjs!(w,"https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.js")
loadjs!(w, katex)
#JQUERY, POPPER, BOOSTRAP JS
customjs = open(f->read(f, String), path*"/ui/scripts/custom.js")
customjs = replace(customjs, "JQUERY"=>path*"/ui/scripts/jquery-3.4.1.slim.min.js")
open(path*"/ui/scripts/custom_tmp.js", "w") do f write(f, customjs) end
loadjs!(w, customjstmp)
loadjs!(w, popper)
loadjs!(w, bsjs)
#ICONS
icons = open(f->read(f, String), path*"/ui/css/icons.css")
icons = replace(icons, "fontseot"=>fontseot)
icons = replace(icons, "fontsttf"=>fontsttf)
icons = replace(icons, "fontswoff"=>fontswoff)
icons = replace(icons, "fontssvg"=>fontssvg)
open(path*"/ui/css/icons_tmp.css", "w") do f write(f, icons) end
loadcss!(w, iconstmp)
# LOAD
index = replace(index, "PHANTOM"=>imphantom)
                     , "SCANNER"=>imscanner)
                     , "PULSES"=>impulses)
body!(w,*(navbar,index,footer))
## UPDATE FUNCTIONS
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
handle(w, "docs") do args...
    content!(w, "#content",open(f->read(f, String), replace(path,"src"=>"")*"docs/build/index.html"))
end
handle(w, "close") do args...
    close(w)
end
## Defaults
global phantom = brain_phantom2D(;axis="coronal")
global seq = []
global scanner = []

w
end

function SpinLab(;frame=false)

## ASSETS
logo = AssetRegistry.register("./src/ui/assets/Logo.png")
loading = AssetRegistry.register("./src/ui/assets/Loading.gif")
jquery = AssetRegistry.register("./src/ui/scripts/jquery-3.4.1.slim.min.js")
popper = AssetRegistry.register("./src/ui/scripts/popper.min.js")
bsjs = AssetRegistry.register("./src/ui/scripts/bootstrap.min.js")
bscss = AssetRegistry.register("./src/ui/css/bootstrap.min.css")
custom = AssetRegistry.register("./src/ui/css/custom.css")
# custom icons
icons = AssetRegistry.register("./src/ui/css/icons.css")
iconstmp = AssetRegistry.register("./src/ui/css/icons_tmp.css")
fontseot = AssetRegistry.register("./src/ui/css/fonts/icomoon.eot")
fontsttf = AssetRegistry.register("./src/ui/css/fonts/icomoon.ttf")
fontswoff = AssetRegistry.register("./src/ui/css/fonts/icomoon.woff")
fontssvg = AssetRegistry.register("./src/ui/css/fonts/icomoon.svg")
## WINDOW
global w = Blink.Window(Dict(
    "title"=>"SpinLab",
    # "x"=>3200,
    # "y"=>150,
    "darkTheme"=>true,
    "autoHideMenuBar"=>true,
    "frame"=>frame, #removes title bar
    "transparent"=>true,
    "backgroundColor"=>"#000",
    "node-integration" => false,
    :icon=>"src/ui/assets/Logo_icon.png",
    "minHeight"=>500,
    "minWidth"=>800
    ))
## LOADING BAR
loadbar = """<center><img src="$loading" width="15%" class="align-middle"></center><br>"""
## NAV BAR
navbar = open(f->read(f, String), "src/ui/html/navbar.html")
navbar = replace(navbar, "LOGO"=>logo)
## CONTENT
index = open(f->read(f, String), "src/ui/html/index.html")
## FOOTER
footer = open(f->read(f, String), "src/ui/html/footer.html")
## Update content
loadcss!(w, bscss)
loadcss!(w, custom)
loadjs!(w, jquery)
loadjs!(w, popper)
icons = open(f->read(f, String), "src/ui/css/icons.css")
icons = replace(icons, "fontseot"=>fontseot)
icons = replace(icons, "fontsttf"=>fontsttf)
icons = replace(icons, "fontswoff"=>fontswoff)
icons = replace(icons, "fontssvg"=>fontssvg)
open("src/ui/css/icons_tmp.css", "w") do f write(f, icons) end
loadcss!(w, iconstmp)
loadjs!(w, bsjs) # boostrap = "<script src=\"$bsjs\"></script>"
body!(w,*(navbar,index,footer))
## UPDATE FUNCTIONS
handle(w, "index") do args...
     @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
     body!(w,*(navbar,index,footer))
end
handle(w, "pulses") do args...
    @js_ w (@var loading = $loadbar; document.getElementById("content").innerHTML=loading)
    include("src/ui/Pulses.jl")
end
handle(w, "close") do args...
    close(w)
end

end

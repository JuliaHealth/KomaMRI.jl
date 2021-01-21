path = @__DIR__

loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Phantom)
map!(f->begin
        print("Loading... $f\n")
        JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"phantom")
    end
    , data, loadbutton)


# Ploting recon
global recon = ifftc(kdata)

hh, ww = 420,550
# l = PlotlyJS.Layout(;title="Reconstruction", yaxis_title="y", yaxis=attr(scaleanchor="x"),
#     xaxis_title="x",height=hh,width=ww,
#     modebar=attr(orientation="v"),scene=attr(aspectratio=attr(x=1,y=1,z=1)))
l = PlotlyJS.Layout(;title="Reconstruction",yaxis_title="y",
    xaxis_title="x",height=hh,width=ww,
    yaxis=attr(scaleanchor="x"),
    modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
p = PlotlyJS.plot(PlotlyJS.heatmap(z=abs.(recon),showscale=false,colorscale="Greys"),l)
plt = Observable{Any}(p)
# PlotlyJS.savefig(p, path*"/assets/phantom.png", width=320, height=300)

function makebuttons(ph)
    global phantom = ph #rewriting default phantom
    prop = propertynames(ph)[4:end-2]
    propnmtuple = string.(prop)
    propnm = [i for i in propnmtuple]
    buttons = button.(propnm)
    for (btn, key, keyname) in zip(buttons, prop, propnm)
        map!(t -> begin
            l = PlotlyJS.Layout(;title="Reconstruction",yaxis_title="y",
                                xaxis_title="x",height=hh,width=ww,
                                yaxis=attr(scaleanchor="x"),
                                modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
            h = PlotlyJS.heatmap(x=ph.x*1e2,y=ph.y*1e2,z=getproperty(ph,key))
            p = PlotlyJS.plot(h,l)
            end
            , plt, btn)
    end
    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)

pulseseq = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", pulseseq)

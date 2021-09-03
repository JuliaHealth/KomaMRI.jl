path = @__DIR__

loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Phantom)
#Assigning function of data when load button (filepicker) is changed
map!(f->begin
        global phantom = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"phantom")
        print("Loaded! ... $f\n")
    end
    , data, loadbutton)

hh, ww = 420,550
# l = PlotlyJS.Layout(;title=phantom.name*": ρ", yaxis_title="y [cm]",
#     yaxis=attr(scaleanchor="x"),
#     xaxis_title="x [cm]",height=hh,width=ww,
#     modebar=attr(orientation="v"),xaxis=attr(constrain="domain"))
# p = PlotlyJS.plot(PlotlyJS.heatmap(x=phantom.x*1e2,y=phantom.y*1e2,
#     z=phantom.ρ),l)
p = @manipulate for t0 = range(0,dur(seq),length=10)*1e3
    l = PlotlyJS.Layout(;title=phantom.name*": ρ",yaxis_title="y [cm]",
        xaxis_title="x [cm]",height=hh,width=ww,
        yaxis=attr(scaleanchor="x"),
        modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
    h = PlotlyJS.scattergl(x=(phantom.x .+ phantom.ux(phantom.x,phantom.y,0,t0*1e-3))*1e2,
                            y=(phantom.y .+ phantom.uy(phantom.x,phantom.y,0, t0*1e-3))*1e2,
                            mode="markers",
                            marker=attr(color=phantom.ρ, showscale=true, colorscale="Viridis"),
                            text=phantom.ρ; 
                            colorbar=attr(ticksuffix=""))
    PlotlyJS.plot(h,l)
end
plt = Observable{Any}(p)

#TODO: Improve this using https://github.com/JuliaGizmos/Interact.jl
function makebuttons(ph)
    global phantom = ph #rewriting default phantom
    prop = propertynames(ph)[5:end-3]
    propnmtuple = string.(prop)
    propnm = [i for i in propnmtuple]
    buttons = button.(propnm)

    for (btn, key, keyname) in zip(buttons, prop, propnm)
        map!(t -> begin
            @manipulate for t0 = range(0,dur(seq),length=10)*1e3
                l = PlotlyJS.Layout(;title=ph.name*": "*keyname,yaxis_title="y [cm]",
                    xaxis_title="x [cm]",height=hh,width=ww,
                    yaxis=attr(scaleanchor="x"),
                    modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
                h = PlotlyJS.scattergl(x=(ph.x .+ ph.ux(ph.x,ph.y,0,t0*1e-3))*1e2,
                                     y=(ph.y .+ ph.uy(ph.x,ph.y,0,t0*1e-3))*1e2,
                                     mode="markers",
                                     marker=attr(color=getproperty(ph,key), showscale=true, colorscale="Viridis"),
                                     text=getproperty(ph,key); 
                                     colorbar=attr(ticksuffix=""))
                p = PlotlyJS.plot(h,l)
            end
        end
        , plt, btn)
    end

    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)
columnbuttons = makebuttons(phantom)
pulseseq = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", pulseseq)

path = @__DIR__

loadbutton = filepicker("Choose .phantom file.."; accept=".phantom", value="")
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Phantom)
hh, ww = 420,550
p = @manipulate for t0 = range(0,dur(seq),length=10)*1e3
    l = PlotlyJS.Layout(;title=phantom.name*": ρ",yaxis_title="y [cm]",
        xaxis_title="x [cm]",height=hh,width=ww,
        yaxis=attr(scaleanchor="x"),
        modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
    h = PlotlyJS.PlotlyJS.scatter3d(x=(phantom.x .+ phantom.ux(phantom.x,phantom.y,phantom.z,t0*1e-3))*1e2,
                                    y=(phantom.y .+ phantom.uy(phantom.x,phantom.y,phantom.z,t0*1e-3))*1e2,
                                    z=(phantom.z .+ phantom.uz(phantom.x,phantom.y,phantom.z,t0*1e-3))*1e2,
                            mode="markers",
                            marker=attr(color=phantom.ρ, showscale=true, colorscale="Viridis"),
                            text=phantom.ρ; 
                            colorbar=attr(ticksuffix=""))
    PlotlyJS.plot(h,l)
end
plt = Observable{Any}(p)
#Assigning function of data when load button (filepicker) is changed
map!(f->begin
        if f!=""
        global pha_file = basename(f)
        global phantom = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"phantom")
        # print("Loaded! ... $f\n")
        else
        phantom
        end
    end
    , data, loadbutton)

function makebuttons(ph)
    global phantom = ph #rewriting default phantom
    prop = propertynames(ph)[5:end-3]
    propnmtuple = string.(prop)
    propnm = [i for i in propnmtuple]
    buttons = button.(propnm)

    for (btn, key, keyname) in zip(reverse(buttons), reverse(prop), reverse(propnm))
        map!(t -> begin
            @manipulate for t0 = range(0,dur(seq),length=10)*1e3
                l = PlotlyJS.Layout(;title=ph.name*": "*keyname,yaxis_title="y [cm]",
                    xaxis_title="x [cm]",height=hh,width=ww,
                    yaxis=attr(scaleanchor="x"),
                    modebar=attr(orientation="v"),xaxis=attr(constrain="domain"),hovermode="closest")
                h = PlotlyJS.scatter3d(x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
                                       y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
                                       z=(ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
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
# columnbuttons = makebuttons(data)
ui = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", ui)

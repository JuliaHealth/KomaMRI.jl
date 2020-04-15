path = @__DIR__

loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(DataFrame)

map!(CSV.read, data, loadbutton)

hh, ww = 420,550
l = PlotlyJS.Layout(;yaxis_title="y [cm]",
    xaxis_title="x [cm]",title=phantom.name*" (p.d.)",height=hh,width=ww,
    modebar=attr(orientation="v"),xaxis=attr(constrain="domain"))
p = PlotlyJS.plot(PlotlyJS.heatmap(x=phantom.x*1e2,y=phantom.y*1e2,
    z=phantom.ρ; colorbar=attr(ticksuffix=" (a.u)")),l)
PlotlyJS.savefig(p, path*"/assets/phantom.png")

plt = Observable{Any}(p)

function makebuttons(df)
    buttons = button.(string.(names(df)))
    for (btn, name) in zip(buttons, names(df))
        map!(t -> begin
            p = PlotlyJS.plot(df[name][:],l)
            p
            end
            , plt, btn)
    end
    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)
pulseseq = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", pulseseq)

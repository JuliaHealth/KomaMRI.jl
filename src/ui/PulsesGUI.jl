loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(DataFrame)

map!(CSV.read, data, loadbutton)

p = PlotlyJS.plot()
hh = 300
PlotlyJS.relayout!(p,height=hh,modebar=attr(orientation="v"))
plt = Observable{Any}(p)

function makebuttons(df)
    buttons = button.(string.(names(df)))
    for (btn, name) in zip(buttons, names(df))
        map!(t -> begin
            p = PlotlyJS.plot(df[name][:])
            PlotlyJS.relayout!(p,height=hh)
            p
            end
            , plt, btn)
    end
    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)
pulseseq = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", pulseseq)

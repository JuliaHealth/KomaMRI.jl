loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Array{Sequence})
map!(f->begin
        print("Loading... $f\n")
        JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"seq")
    end
    , data, loadbutton)

hh, ww = 800, 800
l = PlotlyJS.Layout(;title="Acquired signal", yaxis_title="ky [m^-1]",
    xaxis_title="kx [m^-1]",#height=hh,width=ww,
    modebar=attr(orientation="v"),legend=false)
absS = PlotlyJS.scatter(y=abs.(signal),label="|S(t)|")
reS = PlotlyJS.scatter(y=real.(signal),label="Re{S(t)}")
imS = PlotlyJS.scatter(y=imag.(signal),label="Im{S(t)}")
p = PlotlyJS.plot([absS,reS,imS],l)
plt = Observable{Any}(p)

function makebuttons(df)
    namesdf = 1:length(df)
    buttons = button.(string.())
    for (btn, name) in zip(buttons, namesdf)
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

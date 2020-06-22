loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Array{Sequence})
map!(f->begin
        print("Loading... $f\n")
        JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"seq")
    end
    , data, loadbutton)

hh, ww = 450, 500
l = PlotlyJS.Layout(;title="k-space", yaxis_title="ky [m^-1]",
    xaxis_title="kx [m^-1]",#height=hh,width=ww,
    modebar=attr(orientation="v"),legend=false,height=hh,width=ww)
p = MRIsim.plot_grads(seq)
plt = Observable{Any}(p)

function makebuttons(seq)
    namesseq = ["Sequence","k-space"]
    buttons = button.(string.(namesseq))
    for (btn, name) in zip(buttons, namesseq)
        if name == "Sequence"
            map!(t-> begin
                    MRIsim.plot_grads(seq)
                end
                , plt, btn)
        elseif name == "k-space"
            map!(t->begin
                seqADC = sum([is_DAC_on(s) ? s : Sequence() for s in seq])
                kspace = MRIsim.get_designed_kspace(seqADC)
                line = PlotlyJS.scatter(x=kspace[:,1],y=kspace[:,2],mode="lines")
                marker = PlotlyJS.scatter(x=kspace[:,1],y=kspace[:,2],mode="markers",
                    marker=attr(size=5,color=1:length(kspace),colorscale="Jet"))
                PlotlyJS.plot([line,marker],l)
            end
            , plt, btn)
        end
    end
    dom"div"(hbox(buttons))
end
map!(makebuttons, columnbuttons, data)

columnbuttons = makebuttons(seq)
pulseseq = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", pulseseq)

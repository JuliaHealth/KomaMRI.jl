loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Array{Sequence})
map!(f->begin
        print("Loading... $f\n")
        JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"seq")
    end
    , data, loadbutton)

hh, ww = 600, 600
l = PlotlyJS.Layout(;title=L"k-space", yaxis_title="ky [m^-1]",
    xaxis_title="kx [m^-1]",
    modebar=attr(orientation="v"),height=hh,width=ww,hovermode="closest")
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
                ACQ = is_DAC_on.(seq)
                #Adding gradients before and after acq, to show the centered k-space
                ACQ = (circshift(ACQ,1)+circshift(ACQ,-1)+ACQ).!=0 
                seqADC = sum(seq[ACQ])
                kspace = MRIsim.get_designed_kspace(seqADC)
                N = size(kspace,1)-1

                c = ["hsv(200,$((i-1)/N*255),70)" for i=1:N-1]
                lines = PlotlyJS.scattergl(x=kspace[:,1],y=kspace[:,2],mode="lines+markers",
                        line=attr(color=c),aspectratio=1)
                p = PlotlyJS.plot(lines,l)
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

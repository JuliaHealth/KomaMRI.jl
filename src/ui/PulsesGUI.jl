loadbutton = filepicker("Choose .seq file..."; accept=".seq", value="")
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(Sequence)
dict = Observable{Any}(simParams)
p = plot_seq(seq)
plt = Observable{Any}(p)
#Assigning function of data when load button (filepicker) is changed
map!(f->begin
            if f!=""
                global seq_file = basename(f)
                global seq = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"seq")
            else
            seq
            end
        end
    , data, loadbutton)

function makebuttons(seq)
    namesseq = ["Sequence","k-space"]
    buttons = button.(string.(namesseq))

    for (btn, name) in zip(reverse(buttons), reverse(namesseq))
        if name == "Sequence"
            map!(t-> begin
                plot_seq(seq)
                end
                , plt, btn)
        elseif name == "k-space"
            map!(t->begin
                plot_kspace(seq, simParams)            
            end
            , plt, btn)
        end
    end
    
    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)

# columnbuttons = makebuttons(seq)
ui = dom"div"(loadbutton, columnbuttons, plt)
content!(w, "div#content", ui)
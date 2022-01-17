path = @__DIR__
# loadbutton = filepicker()
# columnbuttons = Observable{Any}(dom"div"())
# data = Observable{Any}(Phantom)
# map!(f->begin
#         global recon = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(f),"recon")
#         print("Loaded! ... $f\n")
#     end
#     , data, loadbutton)
columnbuttons = Observable{Any}(dom"div"())
data = Observable{Any}(signal)
p = plot_image(image)
plt = Observable{Any}(p)
# Ploting recon
function makebuttons(signal)
    namesseq = ["Image","k-data"]
    buttons = button.(string.(namesseq))

    for (btn, name) in zip(reverse(buttons), reverse(namesseq))
        if name == "Image"
            map!(t-> begin
                    recParams["recon"] = "fft"
                    image = reconstruction(signal, recParams)
                    plot_image(image)
                end
                , plt, btn)
        elseif name == "k-data"
            map!(t-> begin
                    recParams["recon"] = "skip"
                    image = reconstruction(signal, recParams)
                    plot_image(image)         
            end
            , plt, btn)
        end
    end
    dom"div"(hbox(buttons))
end

map!(makebuttons, columnbuttons, data)

ui = dom"div"(columnbuttons, plt)
content!(w, "div#content", ui)
#plt = Observable{Any}()
#ui = dom"div"(plt)
## Ploting recon
#map!(t-> begin
#        @manipulate for slice = 1:size(image,3)
#            aux = abs.(image) * prod(size(image)[1:2])
#            plot_image(aux[:,:,slice],zmin=minimum(aux[:]),zmax=maximum(aux[:]);darkmode,title="Reconstruction ($slice/$(size(image,3)))")
#        end
#    end
#    , plt, img_obs)
#content!(w, "div#content", ui)

plt = Observable{Any}()
btn = button("Export .mat")
map!(t-> begin
        @manipulate for slice = 1:size(image,3)
            aux = abs.(image) * prod(size(image)[1:2])
            plot_image(aux[:,:,slice],zmin=minimum(aux[:]),zmax=maximum(aux[:]);darkmode,title="Reconstruction ($slice/$(size(image,3)))")
        end
    end
    , plt, img_obs)
ui = dom"div"(vbox(dom"div"(hbox(btn)), plt))
content!(w, "div#content", ui)

function export2mat()

    matwrite("image.mat", Dict("image" => image))

end

on(n -> export2mat(), btn)

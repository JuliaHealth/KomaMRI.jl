
function plotInteract()
    @manipulate for slice = 1:size(image,3)
        aux = abs.(image) * prod(size(image)[1:2])
        plot_image(aux[:,:,slice],zmin=minimum(aux[:]),zmax=maximum(aux[:]);darkmode,title="Reconstruction ($slice/$(size(image,3)))")
    end
end

plt = Observable{Any}(plotInteract())
map!(t-> plotInteract(), plt, img_obs)
ui = dom"div"(plt)
content!(w, "div#content", ui)
@js_ w document.getElementById("content").dataset.content = "absi"

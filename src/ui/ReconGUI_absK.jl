function plotInteract()
    @manipulate for slice = 1:size(image,3)
        aux = log.(abs.(kspace[:,:,slice]).+1)
        plot_image(aux,zmin=0,zmax=.1*maximum(aux[:]);darkmode,title="Reconstruction ($slice/$(size(image,3)))")
    end
end

plt = Observable{Any}(plotInteract())
map!(t-> plotInteract(), plt, img_obs)
ui = dom"div"(plt)
@js_ w document.getElementById("content").dataset.content = $CONT_ABSK
content!(w, "div#content", ui)

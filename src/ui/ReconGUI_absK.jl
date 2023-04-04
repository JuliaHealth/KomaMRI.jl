plt = Observable{Any}()
map!(t-> begin
        @manipulate for slice = 1:size(image,3)
            aux = log.(abs.(kspace[:,:,slice]).+1)
            plot_image(aux,zmin=0,zmax=.1*maximum(aux[:]);darkmode,title="Reconstruction ($slice/$(size(image,3)))")
        end
    end
    , plt, img_obs)
ui = dom"div"(plt)
content!(w, "div#content", ui)

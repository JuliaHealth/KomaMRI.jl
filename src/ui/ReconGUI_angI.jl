plt = Observable{Any}()
map!(t-> begin
        @manipulate for slice = 1:size(image,3)
            aux = abs.(image) * prod(size(image)[1:2])
            plot_image(angle.(image[:,:,slice]),zmin=-π,zmax=π;darkmode,title="Reconstruction ($slice/$(size(image,3)))")
        end
    end
    , plt, img_obs)
ui = dom"div"(plt)
content!(w, "div#content", ui)

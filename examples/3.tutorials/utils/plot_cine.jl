using Plots, Printf

function plot_cine(frames, fps; Δt=1/fps, filename="cine_recon.gif")

	x = 0:size(frames[1])[2]-1
	y = 1:size(frames[1])[1]

	global_min = minimum(reduce(vcat, frames))  # Mínimo global en todos los frames
    global_max = maximum(reduce(vcat, frames))  # Máximo global en todos los frames

	t = 0

	anim = @animate for image in frames
		t += Δt
		Plots.plot!(
			Plots.heatmap(
				x,y,image',color=:greys; 
				aspect_ratio=:equal, 
				colorbar=true, 
				clim=(global_min, global_max)
			),
			title="t = "*Printf.@sprintf("%.3f", t)*"s", 
			xlims=(minimum(x), maximum(x)), 
			ylims=(minimum(y), maximum(y))
		)
	end

	gif(anim, filename, fps = fps)
end
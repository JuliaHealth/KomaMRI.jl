"""
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)

Define colors for dark or light mode.

# Arguments
- `darkmode`: (`::Bool`) boolean that selects dark or light mode

# Returns
- `bgcolor`: (`::String`) backgound color
- `text_color`: (`::String`) text color
- `plot_bgcolor`: (`::String`) color background for the plots
- `grid_color`: (`::String`) color of the grids
- `sep_color`: (`::String`) color of separator lines
"""
function theme_chooser(darkmode)
	if darkmode
		bgcolor = "rgba(0,0,0,0)"#"rgb(13,16,17)"
		text_color = "gray"
		plot_bgcolor = "rgb(22,26,29)" #"rgb(33,37,41)"
		grid_color = "rgb(40,52,66)" #rgb(40,40,40)
		sep_color = "white"
	else
		bgcolor = "rgba(0,0,0,0)"#"white"
		text_color = "gray"#"rgb(49,70,101)"
		plot_bgcolor = "rgb(229,236,246)"
		grid_color = "white"
		sep_color = "black"
	end
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color
end

"""
    c_map_interp = interp_map(c_map, t_interp)

Interpolates a color map. This is used for plotting the kspace (refer to
[`plot_kspace`](@ref)).

# Arguments
- `c_map`: (`::Vector{Vector{Any}}`)color map. Every element of this vector has a
    vector with a number between 0-1 in its first element and a color string in its second
    element. It serves as a reference to create a color map with more elements
- `t_interp`: (`::Vector{Float64}`) the vector with values between 0-1 that are the
    reference for interpolate the color map with more elements

# Returns
- `c_map_interp`: (`::Vector{String}`) vector with color strings with interpolated
    values
"""
function interp_map(c_map, t_interp)
	idx = [c[1] for c = c_map]
	R = [parse.(Int, split(c[2][5:end-1],", "))[1] for c = c_map]
	G = [parse.(Int, split(c[2][5:end-1],", "))[2] for c = c_map]
	B = [parse.(Int, split(c[2][5:end-1],", "))[3] for c = c_map]
	r = LinearInterpolation(idx,R)(t_interp)
	g = LinearInterpolation(idx,G)(t_interp)
	b = LinearInterpolation(idx,B)(t_interp)
	c_map_interp = ["hsv($r, $g, $b)" for (r,g,b)=zip(r,g,b)]
	c_map_interp
end

"""
    p = plot_seq(seq; width, height, slider, show_seq_blocks, show_sim_blocks, Nblocks,
            darkmode, max_rf_samples, range)

Plots a sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Int64`, `=nothing`) width of the plot
- `height`: (`::Int64`, `=nothing`) height of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to show sequence blocks
- `show_sim_blocks`: (`::Bool`, `=false`) boolean to show simulation blocks
- `Nblocks`: (`::Int64`, `=0`) number of simulation blocks to display
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `max_rf_samples`: (`::Int64`, `=100`) maximum number of RF samples
- `range`: (`::Vector{Float64}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the Sequence struct

# Examples
```julia-repl
julia> seq = read_seq("examples/1.sequences/spiral.seq")

julia> plot_seq(seq)
```
"""
plot_seq(seq::Sequence; width=nothing, height=nothing, slider=true, show_seq_blocks=false, show_sim_blocks=false, Nblocks=0, darkmode=false, max_rf_samples=100, range=[]) = begin
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	idx = ["Gx" "Gy" "Gz"]
	N = length(seq)
	O = size(seq.RF,1)
	ΔT = durs(seq)
	T0 = cumsum([0; ΔT],dims=1)
	off_val = Inf #This removes the unnecessary points in the plot
	#GRADS
	t1x = vcat([get_theo_t(seq.GR[1,i]) .+ T0[i] for i=1:N]...)
	t1y = vcat([get_theo_t(seq.GR[2,i]) .+ T0[i] for i=1:N]...)
	t1z = vcat([get_theo_t(seq.GR[3,i]) .+ T0[i] for i=1:N]...)
	Gx =  vcat([get_theo_A(seq.GR[1,i];off_val) for i=1:N]...)
	Gy =  vcat([get_theo_A(seq.GR[2,i];off_val) for i=1:N]...)
	Gz =  vcat([get_theo_A(seq.GR[3,i];off_val) for i=1:N]...)
	#RFS
	t2 =  vcat([get_theo_t(seq.RF[1,i];max_rf_samples) .+ T0[i] for i=1:N]...)
	R =   vcat([get_theo_A(r;off_val,max_rf_samples) for r = seq.RF]...)
	#ADC
	t3 =  vcat([get_theo_t(seq.ADC[i])  .+ T0[i] for i=1:N]...)
	D =   vcat([get_theo_A(d;off_val) for d = seq.ADC]...)
	#Shapes
	shapes = []
	if show_seq_blocks
		aux = [line(
			xref="x", yref="paper",
			x0=T0[i]*1e3, y0=0,
			x1=T0[i]*1e3, y1=1,
			line=attr(color=sep_color, width=2),
			) for i = 1:N+1]
		append!(shapes, aux)
	end
	# Visually check the simulation blocks
	if show_sim_blocks
		#This is the preparation of the default simulate function
		t, _ = KomaMRI.get_uniform_times(seq, 1e-3)
		breaks = KomaMRI.get_breaks_in_RF_key_points(seq,t)
		Nt = length(t)
		if Nblocks == 0
			Nblocks = ceil(Int, 6506*Nt/1.15e6)
		end
		parts = KomaMRI.kfoldperm(Nt,Nblocks;type="ordered",breaks)
		t_sim_parts = [t[p[1]] for p in parts]
		#Create lines
		aux = [line(
			xref="x", yref="paper",
			x0=t_sim_parts[i]*1e3, y0=0,
			x1=t_sim_parts[i]*1e3, y1=1,
			line=attr(color="Red", width=1),
			) for i = 1:length(t_sim_parts)]
		append!(shapes, aux)
	end
	l = PlotlyJS.Layout(; hovermode="closest",
			xaxis_title="",
			modebar=attr(orientation="h",yanchor="bottom",xanchor="right",y=1,x=0,bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
			legend=attr(orientation="h",yanchor="bottom",xanchor="left",y=1,x=0),
			plot_bgcolor=plot_bgcolor,
			paper_bgcolor=bgcolor,
			xaxis_gridcolor=grid_color,
			yaxis_gridcolor=grid_color,
			xaxis_zerolinecolor=grid_color,
			yaxis_zerolinecolor=grid_color,
			font_color=text_color,
			yaxis_fixedrange = false,
			xaxis=attr(
				ticksuffix=" ms",
				range=range[:],
				rangeslider=attr(visible=slider),
				rangeselector=attr(
					buttons=[
						attr(count=1,
						label="1m",
						step=10,
						stepmode="backward"),
						attr(step="all")
						]
					),
				),
			shapes = shapes,
			margin=attr(t=0,l=0,r=0,b=0)
			)
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	plotter = PlotlyJS.scatter #using scattergl speeds up the plotting but does not show the sequence in the slider below
	p = [plotter() for j=1:(3+2O+1)]
	#GR
	p[1] = plotter(x=t1x*1e3, y=Gx*1e3,name=idx[1],hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
	p[2] = plotter(x=t1y*1e3, y=Gy*1e3,name=idx[2],hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
	p[3] = plotter(x=t1z*1e3, y=Gz*1e3,name=idx[3],hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
	#RF
	for j=1:O
		phase = angle.(R[:,j])
		phase[R[:,j] .== Inf] .= Inf
		p[2j-1+3] = plotter(x=t2*1e3, y=abs.(R[:,j])*1e6,name="|RF_$j|",hovertemplate="(%{x:.4f} ms, %{y:.2f} μT)")
		p[2j+3] =   plotter(x=t2*1e3, y=phase, text=ones(size(t2)), name="<RF_$j",hovertemplate="(%{x:.4f} ms, ∠B1: %{y:.4f} rad)", visible="legendonly")
	end
	#ADC
	p[2O+3+1] = plotter(x=t3*1e3, y=D*1., name="ADC",hovertemplate="(%{x:.4f} ms, %{y:i})")
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)
	PlotlyJS.plot(p, l; config) #, options=Dict(:displayModeBar => false))
end

"""
    p = plot_image(image; height, width, zmin, zmax, darkmode, title)

Plots an image matrix.

# Arguments
- `image`: (`::Matrix{Float64}`) image matrix

# Keywords
- `height`: (`::Int64`, `=750`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `zmin`: (`::Float64`, `=minimum(abs.(image[:]))`) reference value for minimum color
- `zmax`: (`::Float64`, `=maximum(abs.(image[:]))`) reference value for maximum color
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `title`: (`::String`, `=""`) title of the plot

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the image matrix
"""
function plot_image(image; height=600, width=nothing, zmin=minimum(abs.(image[:])), zmax=maximum(abs.(image[:])), darkmode=false, title="")
	#Layout
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = PlotlyJS.Layout(;title=title,yaxis_title="y",
    xaxis_title="x",margin=attr(t=50,l=0,r=0,b=0),
    yaxis=attr(scaleanchor="x"),
	font_color=text_color,
    modebar=attr(orientation="v",bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),xaxis=attr(constrain="domain"),hovermode="closest",
	paper_bgcolor=bgcolor,
	plot_bgcolor=plot_bgcolor,
	xaxis_gridcolor=grid_color,
	yaxis_gridcolor=grid_color,
	xaxis_zerolinecolor=grid_color,
	yaxis_zerolinecolor=grid_color)
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	#Plot
	p = PlotlyJS.heatmap(z=image,transpose=false,zmin=zmin,zmax=zmax,colorscale="Greys")
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)
	PlotlyJS.plot(p,l;config)
end

"""
    p = plot_kspace(seq; width=nothing, height=nothing, darkmode=false)

Plots the k-space of a sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Int64`, `=nothing`) width of the plot
- `height`: (`::Int64`, `=nothing`) height of the plot
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the k-space of the sequence struct `seq`

# Examples
```julia-repl
julia> seq = read_seq("examples/1.sequences/spiral.seq")

julia> plot_kspace(seq)
```
"""
function plot_kspace(seq; width=nothing, height=nothing, darkmode=false)
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	#Calculations of theoretical k-space
	kspace, kspace_adc = get_kspace(seq; Δt=1) #simParams["Δt"])
	t_adc = get_sample_times(seq)
	#Colormap
	c_map = [[t, "hsv($(floor(Int,(1-t)*255)), 100, 50)"] for t=range(0,1;length=10)] # range(s,b,N) only works in Julia 1.7.3
	c = "gray"
	c2_idx = []
	counter = 0
	for s in seq
		if is_ADC_on(s)
			N = s.ADC.N[1]
			append!(c2_idx, counter:N+counter-1)
			counter += N
		end
	end
	c2 = interp_map(c_map, c2_idx ./ maximum(c2_idx))
	#Layout
	mink = minimum(kspace_adc,dims=1)
	maxk = maximum(kspace_adc,dims=1)
	dW = maximum(maxk .- mink, dims=2) * .3
	mink .-= dW
	maxk .+= dW
	#Layout
	l = PlotlyJS.Layout(;
		paper_bgcolor=bgcolor,
		scene=attr(xaxis=attr(title="kx [m⁻¹]",range=[mink[1],maxk[1]],backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
				   yaxis=attr(title="ky [m⁻¹]",range=[mink[2],maxk[2]],backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
				   zaxis=attr(title="kz [m⁻¹]",range=[mink[3],maxk[3]],backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color)),
		modebar=attr(orientation="h",yanchor="bottom",xanchor="right",y=1,x=0,bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
		legend=attr(orientation="h",yanchor="bottom",xanchor="left",y=1,x=0),
		font_color=text_color,
		scene_camera_eye=attr(x=0, y=0, z=1.7),
		scene_camera_up=attr(x=0, y=1., z=0),
		scene_aspectmode="cube",
		margin=attr(t=0,l=0,r=0))
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	#Plot
	p = [PlotlyJS.scatter() for j=1:3]
	p[1] = PlotlyJS.scatter3d(x=kspace[:,1],y=kspace[:,2],z=kspace[:,3],mode="lines",
			line=attr(color=c),name="Trajectory",hoverinfo="skip")
	p[2] = PlotlyJS.scatter3d(x=kspace_adc[:,1],y=kspace_adc[:,2],z=kspace_adc[:,3],text=round.(t_adc*1e3,digits=3),mode="markers",
			line=attr(color=c2),marker=attr(size=2),name="ADC",hovertemplate="kx: %{x:.1f} m⁻¹<br>ky: %{y:.1f} m⁻¹<br>kz: %{z:.1f} m⁻¹<br><b>t_acq</b>: %{text} ms<extra></extra>")
	p[3] = PlotlyJS.scatter3d(x=[0],y=[0],z=[0],name="k=0",marker=attr(symbol="cross",size=10,color="red"))
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	PlotlyJS.plot(p,l; config)
end

"""
    p = plot_M0(seq; height=nothing, width=nothing, slider=true, darkmode=false)

Plots the magnetization M0 of a sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `height`: (`::Int64`, `=nothing`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the magnetization M0 of the sequence struct `seq`

# Examples
```julia-repl
julia> seq = read_seq("examples/1.sequences/spiral.seq")

julia> plot_M0(seq)
```
"""
function plot_M0(seq; height=nothing, width=nothing, slider=true, darkmode=false)
	#Times
	dt = 1
	t, Δt = KomaMRI.get_uniform_times(seq, dt)
	#kx,ky
	ts = t .+ Δt
	rf_idx, rf_type = KomaMRI.get_RF_types(seq, t)
	k, _ =  KomaMRI.get_kspace(seq; Δt=dt)

	#plots k(t)
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = PlotlyJS.Layout(;yaxis_title="M0", hovermode="closest",
			xaxis_title="",
			modebar=attr(orientation="h",yanchor="bottom",xanchor="right",y=1,x=0,bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
			legend=attr(orientation="h",yanchor="bottom",xanchor="left",y=1,x=0),
			plot_bgcolor=plot_bgcolor,
			paper_bgcolor=bgcolor,
			xaxis_gridcolor=grid_color,
			yaxis_gridcolor=grid_color,
			xaxis_zerolinecolor=grid_color,
			yaxis_zerolinecolor=grid_color,
			font_color=text_color,
			yaxis_fixedrange = false,
			xaxis=attr(
				ticksuffix=" ms",
				range=[0.,min(20,1e3*dur(seq))],
				rangeslider=attr(visible=slider),
				rangeselector=attr(
					buttons=[
						attr(count=1,
						label="1m",
						step=10,
						stepmode="backward"),
						attr(step="all")
						]
					),
				),
			margin=attr(t=0,l=0,r=0)
			)
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	plotter = PlotlyJS.scatter #using scattergl speeds up the plotting but does not show the sequence in the slider below
	p = [plotter() for j=1:4]
	p[1] = plotter(x=ts*1e3, y=k[:,1], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)", name="x")
	p[2] = plotter(x=ts*1e3, y=k[:,2], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)", name="y")
	p[3] = plotter(x=ts*1e3, y=k[:,3], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)", name="z")
	p[4] = plotter(x=t[rf_idx]*1e3,y=rf_type,name="RFs",marker=attr(symbol="cross",size=8,color="orange"),mode="markers")
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "select", "select2d","lasso2d", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)
	PlotlyJS.plot(p, l; config)
end

"""
    p = plot_phantom_map(ph, key; t0=0, height=600, width=nothing, darkmode=false)

Plots a phantom map for a specific spin parameter given by `key`.

# Arguments
- `ph`: (`::Phantom`) Phantom struct
- `key`: (`::Symbol`, opts: [`:ρ`, `:T1`, `:T2`, `:T2s`, `:x`, `:y`, `:z`]) symbol for
    displaying different parameters of the phantom spins

# Keywords
- `t0`: (`::Float64`, `=0`, `[ms]`) time to see displacement of the phantom
- `height`: (`::Int64`, `=600`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the phantom map for a specific spin parameter

# References
Colormaps from https://github.com/markgriswold/MRFColormaps
Towards Unified Colormaps for Quantitative MRF Data, Mark Griswold, et al. (2018).

# Examples
```julia-repl
julia> obj2D, obj3D = brain_phantom2D(), brain_phantom3D();

julia> plot_phantom_map(obj2D, :ρ)

julia> plot_phantom_map(obj3D, :ρ)
```
"""
function plot_phantom_map(ph::Phantom, key::Symbol; t0=0, height=600, width=nothing, darkmode=false)
	path = @__DIR__
	cmin_key = minimum(getproperty(ph,key))
	cmax_key = maximum(getproperty(ph,key))
	if key == :T1 || key == :T2 || key == :T2s
		cmin_key = 0
		factor = 1e3
		unit = " ms"
		if key  == :T1
			cmax_key = 2500/factor
			colors = MAT.matread(path*"/assets/T1cm.mat")["T1colormap"]
			N, _ = size(colors)
			idx = range(0,1;length=N) #range(0,T,N) works in Julia 1.7
			colormap = [[idx[n], "rgb($(floor(Int,colors[n,1]*255)),$(floor(Int,colors[n,2]*255)),$(floor(Int,colors[n,3]*255)))"] for n=1:N]
		elseif key == :T2 || key == :T2s
			if key == :T2
				cmax_key = 250/factor
			end
    		colors = MAT.matread(path*"/assets/T2cm.mat")["T2colormap"]
			N, _ = size(colors)
			idx = range(0,1;length=N) #range(0,T,N) works in Julia 1.7
			colormap = [[idx[n], "rgb($(floor(Int,colors[n,1]*255)),$(floor(Int,colors[n,2]*255)),$(floor(Int,colors[n,3]*255)))"] for n=1:N]
		end
	elseif key == :x || key == :y || key == :z
		factor = 1e2
		unit = " cm"
		colormap="Greys"
	elseif key == :Δw
		factor = 1/(2π)
		unit = " Hz"
		colormap="Greys"
	else
		factor = 1
		cmin_key = 0
		unit=""
		colormap="Greys"
	end
	cmin_key *= factor
	cmax_key *= factor
	x0 = minimum([ph.x;ph.y;ph.z])*1e2
    xf = maximum([ph.x;ph.y;ph.z])*1e2
	#Layout
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = PlotlyJS.Layout(;title=ph.name*": "*string(key), margin=attr(t=50,l=0,r=0,b=0),
		paper_bgcolor=bgcolor,font_color=text_color,
		scene=attr(
			xaxis=attr(title="x",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			yaxis=attr(title="y",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			zaxis=attr(title="z",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
		),
		modebar=attr(orientation="h",bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),xaxis=attr(constrain="domain"),
		hovermode="closest",aspectmode="cube")
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	h = PlotlyJS.scatter3d( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							z=(ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							mode="markers",
							marker=attr(color=getproperty(ph,key)*factor,
										showscale=true,
										colorscale=colormap,
										colorbar=attr(ticksuffix=unit, title=string(key)),
										cmin=cmin_key,
										cmax=cmax_key,
										size=2
										),
							text=round.(getproperty(ph,key)*factor,digits=4),
							hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br>z: %{z:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	p = PlotlyJS.plot(h,l;config)
end

"""
    p = plot_signal(raw::RawAcquisitionData; height, width, slider, show_sim_blocks,
            darkmode, range)

Plots a raw signal in ISMRMRD format.

# Arguments
- `raw`: (`::RawAcquisitionData`) RawAcquisitionData struct which is the raw signal in
    ISMRMRD format

# Keywords
- `width`: (`::Int64`, `=nothing`) width of the plot
- `height`: (`::Int64`, `=nothing`) height of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `show_sim_blocks`: (`::Bool`, `=false`) boolean to show simulation blocks
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `range`: (`::Vector{Float64}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the raw signal

# Examples
```julia-repl
julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq("examples/1.sequences/spiral.seq");

julia> ismrmrd = simulate(obj, seq, sys);

julia> plot_signal(ismrmrd)
```
"""
function plot_signal(raw::RawAcquisitionData; width=nothing, height=nothing, slider=true, show_sim_blocks=false, darkmode=false, range=[])
	not_Koma = raw.params["systemVendor"] != "KomaMRI.jl"
	t = []
	signal = []
	current_t0 = 0
	for p in raw.profiles
		dt = p.head.sample_time_us != 0 ? p.head.sample_time_us * 1e-3 : 1
		t0 = p.head.acquisition_time_stamp * 1e-3 #This parameter is used in Koma to store the time offset
		N =  p.head.number_of_samples != 0 ? p.head.number_of_samples : 1
		if not_Koma
			t0 = current_t0 * dt
			current_t0 += N
		end
		append!(t, t0.+(0:dt:dt*(N-1)))
		append!(signal, p.data[:,1]) #Just one coil
		#To generate gap
		append!(t, t[end])
		append!(signal, [Inf + Inf*1im])
	end
	#Show simulation blocks
	shapes = []
	if !not_Koma && show_sim_blocks
		t_sim_parts = raw.params["userParameters"]["t_sim_parts"]
		aux = [line(
			xref="x", yref="paper",
			x0=t_sim_parts[i]*1e3, y0=0,
			x1=t_sim_parts[i]*1e3, y1=1,
			line=attr(color="Red", width=1),
			) for i = 1:length(t_sim_parts)]
		append!(shapes, aux)
	end
	#PLOT
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = PlotlyJS.Layout(; hovermode="closest",
			xaxis_title="",
			modebar=attr(orientation="h",yanchor="bottom",xanchor="right",y=1,x=0,bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
			legend=attr(orientation="h",yanchor="bottom",xanchor="left",y=1,x=0),
			plot_bgcolor=plot_bgcolor,
			paper_bgcolor=bgcolor,
			xaxis_gridcolor=grid_color,
			yaxis_gridcolor=grid_color,
			xaxis_zerolinecolor=grid_color,
			yaxis_zerolinecolor=grid_color,
			font_color=text_color,
			yaxis_fixedrange = false,
			xaxis=attr(
				ticksuffix=" ms",
				range=range[:],
				rangeslider=attr(visible=slider),
				rangeselector=attr(
					buttons=[
						attr(count=1,
						label="1m",
						step=10,
						stepmode="backward"),
						attr(step="all")
						]
					),
				),
			shapes = shapes,
			margin=attr(t=0,l=0,r=0,b=0)
			)
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	plotter = PlotlyJS.scatter
	p = [plotter() for j=1:3]
	p[1] = plotter(x=t,y=abs.(signal), name="|S(t)|",hovertemplate="(%{x:.4f} ms, %{y:.3f} a.u.)")
	p[2] = plotter(x=t,y=real.(signal),name="Re{S(t)}",hovertemplate="(%{x:.4f} ms, %{y:.3f} a.u.)")
	p[3] = plotter(x=t,y=imag.(signal),name="Im{S(t)}",hovertemplate="(%{x:.4f} ms, %{y:.3f} a.u.)")
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
		# modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)
	PlotlyJS.plot(p, l;config)
end

"""
    str = plot_dict(dict::Dict)

Generates a string in html format of the dictionary `dict`.

# Arguments
- `dict`: (`::Dict`) dictionary to generate the html string

# Returns
- `str`: (`::String`) string of the dictionary `dict` which is a table in html format
"""
function plot_dict(dict::Dict)
	html = """
	<table class="table table-dark table-striped">
		<thead>
			<tr>
			<th scope="col">#</th>
			<th scope="col">Name</th>
			<th scope="col">Value</th>
			</tr>
		</thead>
		<tbody>
	"""
	i = 1
	for (key,val) = dict
		html *= """
			<tr>
				<th scope="row">$i</th>
				<td>$(string(key))</td>
				<td>$(string(val))</td>
			</tr>
		"""
		i += 1
	end
	html *= "</tbody></table>"
end




# """Plots gradient moments M0, M1 and M2 for specified axis. Does NOT consider RF pulses."""
# plot_grads_moments(seq::Sequence, idx::Int=1; title="", mode="quick") = begin
# 	if mode == "quick"
# 		plotter = PlotlyJS.scattergl
# 	else
# 		plotter = PlotlyJS.scatter
# 	end
# 	names = ["Gx", "Gy", "Gz"]
# 	M, N = size(seq.GR)
# 	G = [seq.GR[j,floor(Int,i/2)+1].A for i=0:2*N-1, j=1:M]
# 	T = [seq.GR[1,i].T for i=1:N]
# 	t = [sum(T[1:i]) for i=1:N]
# 	t = [t[floor(Int,i/2)+1] for i=0:2*N-1]
# 	t = [0; t[1:end-1]]

# 	M0, M1, M2 = get_M0_M1_M2(seq)
# 	M0t = [M0(t)[idx] for t=t][:]; M0max = maximum(abs.(M0t))*1.1e6
# 	M1t = [M1(t)[idx] for t=t][:]; M1max = maximum(abs.(M1t))*1.1e9
# 	M2t = [M2(t)[idx] for t=t][:]; M2max = maximum(abs.(M2t))*1.1e12
# 	Gmax = maximum(abs.(G))*1.1e3

# 	p = [plotter() for k = 1:4]

# 	l = PlotlyJS.Layout(;
# 	title=title,
# 	xaxis=attr(domain = [0, .75]),
# 	yaxis=attr(title=attr(text="G [mT/m]", standoff=0),
# 		side="left", tickfont=attr(color="#636efa"), titlefont=attr(color="#636efa"),range=[-Gmax, Gmax]),
# 	yaxis2=attr(title=attr(text="M0 [mT/m⋅ms]", standoff=1, position=1), showgrid=false,
# 		anchor="free", overlaying="y",position=0.75, side="right",
# 		tickfont=attr(color="#ef553b"), titlefont=attr(color="#ef553b"), range=[-M0max, M0max]),
# 	yaxis3=attr(title=attr(text="M1 [mT/m⋅ms²]", standoff=1), showgrid=false,
# 		anchor="free",overlaying="y",position=0.75+.25/3,side="right",
# 		tickfont=attr(color="#45d9b2"), titlefont=attr(color="#45d9b2"), range=[-M1max, M1max]),
# 	yaxis4=attr(title=attr(text="M2 [mT/m⋅ms³]", standoff=1), showgrid=false,
# 		anchor="free",overlaying="y",position=0.75+.5/3,side="right",
# 		tickfont=attr(color="#b373fa"), titlefont=attr(color="#b373fa"), range=[-M2max, M2max]),
# 	xaxis_title="t [ms]", height=300)
# 	p[1] = plotter(x=t*1e3, y=G[:,idx]*1e3, name=names[idx] , line_shape="hv")
# 	p[2] = plotter(x=t*1e3, y=M0t*1e6, name="M0", yaxis="y2",line=attr(dash="dash"))
# 	p[3] = plotter(x=t*1e3, y=M1t*1e9, name="M1", yaxis="y3",line=attr(dash="dash"))
# 	p[4] = plotter(x=t*1e3, y=M2t*1e12, name="M2", yaxis="y4",line=attr(dash="dash"))
# 	PlotlyJS.plot(p, l)
# end
# plot_Phantom(obj::Phantom,filename::String) = begin
# 	# Phantom
# 	p1 = heatmap(obj.x*1e2,obj.y*1e2,obj.ρ,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("Proton density")
# 	p2 = heatmap(obj.x*1e2,obj.y*1e2,obj.T2*1e3,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("T2 [ms]")
# 	p3 = heatmap(obj.x*1e2,obj.y*1e2,obj.Δw/2π,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("df [Hz]")
# 	P(i,j) = rotz(Dθ[i,j])[1:2,1:2]; D(i,j) = [obj.Dλ1[i,j] 0;0 obj.Dλ2[i,j]]
# 	nx = [1;0]; ny = [0;1]
# 	Dx = [nx'*P(i,j)'*D(i,j)*P(i,j)*nx for i=1:size(obj.Dλ1,1),j=1:size(obj.Dλ1,2)]
# 	Dy = [ny'*P(i,j)'*D(i,j)*P(i,j)*ny for i=1:size(obj.Dλ1,1),j=1:size(obj.Dλ1,2)]
# 	p4 = heatmap(obj.x*1e2,obj.y*1e2,Dx*1e12,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("Dx [um2/s]")
# 	p5 = heatmap(obj.x*1e2,obj.y*1e2,Dy*1e12,aspect_ratio=:equal)#,clims=(0,maximum(T2[:])))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!("Dy [um2/s]")
# 	p = plot(p1,p2,p3,p4,p5,size=(1300,500),layout=@layout [a b c; d e])
# 	savefig(p,filename)
# end
# plot_sim_res(obj::Phantom,SEQ::Array{Sequence},S::Array{ComplexF64},
# 	t::Array{Float64,2},filename::String,t_k0::Float64,Ga::Float64) = begin
# 	T2m = sum(obj.T2.*obj.ρ)/sum(obj.ρ)
# 	b, n = get_bvalue(SEQ[1])
# 	println("Sequence with b-value: "*string(round(b*1e-6,digits=2))*" s/mm2") # s/mm2
# 	P(i) = rotz(obj.Dθ[i])[1:2,1:2]; D(i) = [obj.Dλ1[i] 0;0 obj.Dλ2[i]]
# 	Deq = [n'*P(i)'*D(i)*P(i)*n for i=1:prod(size(obj.Dλ1))]
# 	Dm = sum(Deq.*obj.ρ)/sum(obj.ρ) #MODIFYY

# 	p = plot_grads(SEQ,t,t_k0,Ga)
# 	p3 = plot([t_k0; t_k0]*1e3,[0; 1.1],linewidth=0.25,color=:black,label="k=0")
# 	plot!(t[:]*1e3,abs.(S),label="|S|")
# 	plot!(t[:]*1e3,exp.(-t[:]/T2m),linestyle=:dash,label="exp(-t/T2)")
# 	# exp(-b*D) <-> exp(-4*π^2*(Δ-δ/3)*q'*D*q)
# 	plot!(t[:]*1e3,exp.(-t[:]/T2m.-b*Dm).*(t[:].>=(Δ+δ)),linestyle=:dash,label="exp(-t/T2-bD)")

# 	xlabel!("Time [ms]")
# 	ylabel!("Signal [a.u]")
# 	xlims!((minimum(t),minimum(t).+dur(sum(SEQ))).*1e3)
# 	p = plot(p[1],p[2],p3,size=(800,600),layout=@layout [a ; b; c])
# 	savefig(p,filename)
# end
# plot_ksapce_trajectory(ACQ::Sequence,t::Array{Float64,2},filename::String) = begin
# 	k = get_designed_kspace(ACQ)
# 	p = plot(k[:,1],k[:,2],legend=:none,aspect_ratio=:equal)
# 	k = get_actual_kspace(ACQ,t)
# 	scatter!(k[:,1],k[:,2],legend=:none,markersize=1)
# 	xlabel!("kx [1/m]"); ylabel!("ky [1/m]")
# 	savefig(p,filename)
# end
# plot_recon(kdata::Array{ComplexF64},rec::Array{ComplexF64},
# 	Δx_pix::Float64,Δy_pix::Float64,filename::String,title::String) = begin
# 	Nx, Ny = size(rec)
# 	xr = -Δx_pix*(Nx-1)/2:Δx_pix:Δx_pix*(Nx-1)/2
# 	yr = -Δy_pix*(Ny-1)/2:Δy_pix:Δy_pix*(Ny-1)/2
# 	Wx, Wy = 1/Δx_pix, 1/Δy_pix
# 	kx = range(-Wx/2,stop=Wx/2,length=Nx)
# 	ky = range(-Wy/2,stop=Wy/2,length=Ny)
# 	p1 = heatmap(kx,ky,log.(abs.(kdata).+1),aspect_ratio=:equal,
# 				legend=:none,size=(400,400))
# 	xaxis!("kx [1/m]"); yaxis!("ky [1/m]")
# 	title!("k-space")
# 	p2 = heatmap(xr*1e2,yr*1e2,abs.(rec),aspect_ratio=:equal,
# 				legend=:none,size=(400,400))
# 	xaxis!("x [cm]"); yaxis!("y [cm]")
# 	title!(title)
# 	#savefig(p1,filename*"_ksp.pdf")
# 	savefig(p2,filename*".pdf")
# end
# plot_Eq(vx,vy,S,S0,SEQ) = begin
# 	Nq, Nθ = size(SEQ)
# 	rs = range(0,stop=1,length=Nq)
# 	θs = range(0,stop=π,length=Nθ); θs = [θs[1:end-1]; π.+θs[1:end-1]];
#     Eq(vx,vy) = begin
#         Eq = [S[i,j][vx,vy] for i=1:Nq,j=1:Nθ]
#         Eq = [Eq[:,1:end-1] Eq[:,1:end-1]]
#     end
#     c = RGB{Float64}(1,1,1)
#     α = real.(abs.(S0[vx,vy]))/maximum(real.(abs.(S0)))
#     pyplot()
#     hm = heatmap(θs,rs,Eq(vx,vy)*(α<.2 ? 0 : 1),proj=:polar,aspect_ratio=:equal,
#     	legend=:none,grid=false,xticks=:none,yticks=:none,
#     	background_color_subplot=α*c)
# end
# plot_Pr(vx,vy,S,S0,SEQ) = begin
# 	Nq, Nθ = size(SEQ)
# 	rs = range(0,stop=1,length=10)*25e-6
# 	θs = range(0,stop=2π,length=10)
#     Pr(vx,vy) = begin
#         Eq = [S[i,j][vx,vy] for i=1:Nq,j=1:Nθ]
# 		b = [get_bvalue(SEQ[i,j])[1] for i=1:Nq,j=1:Nθ]
# 		n = [get_bvalue(SEQ[i,j])[2] for i=1:Nq,j=1:Nθ]
# 		B = [-b[i,j]*[n[i,j][1]^2
# 					2*n[i,j][1]*n[i,j][2]
# 					  n[i,j][2]^2] for i=1:Nq,j=1:Nθ]
# 		B = [B[i][j] for i=1:Nq*Nθ,j=1:3] #reshape
# 		y = log.(abs.(Eq))[:]
# 		# exp(-b*n'D n) <-> exp(-4*π^2*(Δ-δ/3)*q'*D*q)
# 		# min D || -b*n'D*n - log(S) ||_2 + λ ||D||_2
# 		# min D || B D - y ||_2 + λ ||D||_2 => D = (B'B+λI)^-1 B' y
# 		λ = 0
# 		Id = Matrix{Float64}(I,3,3); D = [0;0;0]
# 		for n = 0:2 #Tikhonov regularized and iterativly weighted
# 			W = n==0 ? diagm(0=>Eq[:].^2) : diagm(0=>exp.(2*B*D))
# 			D = (B'*W*B + λ*Id)^-1*B'*W*y
# 		end
# 		# Diffusion propagator
# 		D_inv = [D[1] D[2];
# 				 D[2] D[3]]^-1
# 		pr = exp.(-[([r*cos(θ) r*sin(θ)]*D_inv*[r*cos(θ);r*sin(θ)])[1] for r=rs, θ=θs])
#     end
#     c = RGB{Float64}(1,1,1)
#     α = real.(abs.(S0[vx,vy]))/maximum(real.(abs.(S0)))
#     pyplot()
#     hm = heatmap(θs,rs,Pr(vx,vy)*(α<.2 ? 0 : 1),proj=:polar,aspect_ratio=:equal,
#     	legend=:none,grid=false,xticks=:none,yticks=:none,
#     	background_color_subplot=α*c)
# end

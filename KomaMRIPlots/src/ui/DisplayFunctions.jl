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

function generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0)
	#LAYOUT
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
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
	l = Layout(;title=title, hovermode="closest",
			xaxis_title="",
			modebar=attr(orientation="h", yanchor="bottom", xanchor="right", y=1, x=0, bgcolor=bgcolor, color=text_color, activecolor=plot_bgcolor),
			legend=attr(orientation="h", yanchor="bottom", xanchor="left", y=1, x=0),
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
	#CONFIG
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan", 
								"tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)

	l, config
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
	r = linear_interpolation(idx,R)(t_interp)
	g = linear_interpolation(idx,G)(t_interp)
	b = linear_interpolation(idx,B)(t_interp)
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
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_seq(seq)
```
"""
plot_seq(seq::Sequence; width=nothing, height=nothing, slider=true, show_seq_blocks=false, darkmode=false, max_rf_samples=100, range=[], title="") = begin
	idx = ["Gx" "Gy" "Gz"]
	N = length(seq)
	O = size(seq.RF,1)
	ΔT = KomaMRICore.durs(seq)
	T0 = cumsum([0; ΔT],dims=1)
	off_val = Inf #This removes the unnecessary points in the plot
	#GRADS
	t1x = vcat([KomaMRICore.get_theo_t(seq.GR[1,i]) .+ T0[i] for i=1:N]...)
	t1y = vcat([KomaMRICore.get_theo_t(seq.GR[2,i]) .+ T0[i] for i=1:N]...)
	t1z = vcat([KomaMRICore.get_theo_t(seq.GR[3,i]) .+ T0[i] for i=1:N]...)
	Gx =  vcat([KomaMRICore.get_theo_A(seq.GR[1,i];off_val) for i=1:N]...)
	Gy =  vcat([KomaMRICore.get_theo_A(seq.GR[2,i];off_val) for i=1:N]...)
	Gz =  vcat([KomaMRICore.get_theo_A(seq.GR[3,i];off_val) for i=1:N]...)
	#RFS
	t2 =  vcat([KomaMRICore.get_theo_t(seq.RF[1,i];max_rf_samples) .+ T0[i] for i=1:N]...)
	R =   vcat([KomaMRICore.get_theo_A(r;off_val,max_rf_samples) for r = seq.RF]...)
	#ADC
	t3 =  vcat([KomaMRICore.get_theo_t(seq.ADC[i])  .+ T0[i] for i=1:N]...)
	D =   vcat([KomaMRICore.get_theo_A(d;off_val) for d = seq.ADC]...)
<<<<<<< HEAD
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
		t, _ = KomaMRICore.get_uniform_times(seq, 1e-3)
		breaks = KomaMRICore.get_breaks_in_RF_key_points(seq,t)
		Nt = length(t)
		if Nblocks == 0 #TODO: This should change to a call to a function that generates the default parameters for the simulation
			Nblocks = 20
		end
		parts = KomaMRICore.kfoldperm(Nt,Nblocks;type="ordered",breaks)
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
	l = PlotlyJS.Layout(;title=title, hovermode="closest",
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
=======
    #Plot
	p = [scatter() for j=1:(3+2O+1)]
>>>>>>> koma_v074
	#GR
	p[1] = scatter(x=t1x*1e3, y=Gx*1e3,name=idx[1],hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
	p[2] = scatter(x=t1y*1e3, y=Gy*1e3,name=idx[2],hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
	p[3] = scatter(x=t1z*1e3, y=Gz*1e3,name=idx[3],hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
	#RF
	for j=1:O
		phase = angle.(R[:,j])
		phase[R[:,j] .== Inf] .= Inf
		p[2j-1+3] = scatter(x=t2*1e3, y=abs.(R[:,j])*1e6,name="|RF_$j|",hovertemplate="(%{x:.4f} ms, %{y:.2f} μT)")
		p[2j+3] =   scatter(x=t2*1e3, y=phase, text=ones(size(t2)), name="<RF_$j",hovertemplate="(%{x:.4f} ms, ∠B1: %{y:.4f} rad)", visible="legendonly")
	end
	#ADC
	p[2O+3+1] = scatter(x=t3*1e3, y=D*1., name="ADC",hovertemplate="(%{x:.4f} ms, %{y:i})")
	#Layout and config
	l, config = generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0)
	plot(p, l; config)
end

"""
    p = plot_M0(seq; height=nothing, width=nothing, slider=true, darkmode=false, range=[])

Plots the zero order moment (M0) of a Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `height`: (`::Int64`, `=nothing`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `range`: (`::Vector{Float64}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the moment M0 of the sequence struct `seq`

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_M0(seq)
```
"""
function plot_M0(seq; width=nothing, height=nothing, slider=true, show_seq_blocks=false, darkmode=false, range=[], title="")
	#Times
	dt = 1
	t, Δt = KomaMRICore.get_uniform_times(seq, dt)
	t = t[1:end-1]
	ΔT = KomaMRICore.durs(seq)
	T0 = cumsum([0; ΔT],dims=1)
	#M0
	ts = t .+ Δt
	rf_idx, rf_type = KomaMRICore.get_RF_types(seq, t)
	k, _ =  KomaMRICore.get_kspace(seq; Δt=dt)
	#plots M0
	p = [scatter() for j=1:4]
	p[1] = scatter(x=ts*1e3, y=k[:,1], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)", name="M0x")
	p[2] = scatter(x=ts*1e3, y=k[:,2], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)", name="M0y")
	p[3] = scatter(x=ts*1e3, y=k[:,3], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms)", name="M0z")
	p[4] = scatter(x=t[rf_idx]*1e3,y=rf_type,name="RFs",marker=attr(symbol="cross",size=8,color="orange"),mode="markers")
	#Layout and config
	l, config = generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0)
	plot(p, l; config)
end

"""
    p = plot_M1(seq; height=nothing, width=nothing, slider=true, darkmode=false, range=[])

Plots the first order moment (M1) of a Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence

# Keywords
- `height`: (`::Int64`, `=nothing`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `range`: (`::Vector{Float64}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the moment M1 of the sequence struct `seq`

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_M1(seq)
```
"""
function plot_M1(seq; width=nothing, height=nothing, slider=true, show_seq_blocks=false, darkmode=false, range=[], title="")
	#Times
	dt = 1
	t, Δt = KomaMRICore.get_uniform_times(seq, dt)
	t = t[1:end-1]
	ΔT = KomaMRICore.durs(seq)
	T0 = cumsum([0; ΔT],dims=1)
	#M1
	ts = t .+ Δt
	rf_idx, rf_type = KomaMRICore.get_RF_types(seq, t)
	k, _ =  KomaMRICore.get_M1(seq; Δt=dt)
	#plots M1
	p = [scatter() for j=1:4]
	p[1] = scatter(x=ts*1e3, y=k[:,1], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms²)", name="M1x")
	p[2] = scatter(x=ts*1e3, y=k[:,2], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms²)", name="M1y")
	p[3] = scatter(x=ts*1e3, y=k[:,3], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms²)", name="M1z")
	p[4] = scatter(x=t[rf_idx]*1e3,y=rf_type,name="RFs",marker=attr(symbol="cross",size=8,color="orange"),mode="markers")
	#Layout and config
	l, config = generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0)
	plot(p, l; config)
end


"""
    p = plot_M2(seq; height=nothing, width=nothing, slider=true, darkmode=false, range=[])

Plots the second order moment (M2) of a Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence

# Keywords
- `height`: (`::Int64`, `=nothing`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `range`: (`::Vector{Float64}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the moment M2 of the sequence struct `seq`

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_M2(seq)
```
"""
function plot_M2(seq; width=nothing, height=nothing, slider=true, show_seq_blocks=false, darkmode=false, range=[], title="")
	#Times
	dt = 1
	t, Δt = KomaMRICore.get_uniform_times(seq, dt)
	t = t[1:end-1]
	ΔT = KomaMRICore.durs(seq)
	T0 = cumsum([0; ΔT],dims=1)
	#M2
	ts = t .+ Δt
	rf_idx, rf_type = KomaMRICore.get_RF_types(seq, t)
	k, _ =  KomaMRICore.get_M2(seq; Δt=dt)
	#Plor M2
	p = [scatter() for j=1:4]
	p[1] = scatter(x=ts*1e3, y=k[:,1], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)", name="M2x")
	p[2] = scatter(x=ts*1e3, y=k[:,2], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)", name="M2y")
	p[3] = scatter(x=ts*1e3, y=k[:,3], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)", name="M2z")
	p[4] = scatter(x=t[rf_idx]*1e3,y=rf_type,name="RFs",marker=attr(symbol="cross",size=8,color="orange"),mode="markers")
	#Layout and config
	l, config = generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0)
	plot(p, l; config)
end


"""
    p = plot_eddy_currents(seq, λ; α=ones(size(λ)), height=nothing, width=nothing, slider=true, darkmode=false, range=[])

Plots the eddy currents of a Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence
- `λ`: (`::Float64`, `[s]`) eddy currents decay constant time

# Keywords
- `α`: (`::Vector{Float64}`, `=ones(size(λ))`) eddy currents factors
- `height`: (`::Int64`, `=nothing`) height of the plot
- `width`: (`::Int64`, `=nothing`) width of the plot
- `slider`: (`::Bool`, `=true`) boolean to display a slider
- `darkmode`: (`::Bool`, `=false`) boolean to define colors for darkmode
- `range`: (`::Vector{Float64}`, `=[]`) time range to be displayed initially

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the eddy currents of the sequence struct `seq`

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_eddy_currents(seq, 80e-3)
```
"""
function plot_eddy_currents(seq, λ; α=ones(size(λ)), width=nothing, height=nothing, slider=true, show_seq_blocks=false, darkmode=false, range=[], title="")
	#Times
	dt = 1
	t, Δt = KomaMRICore.get_uniform_times(seq, dt)
	t = t[1:end-1]
	ΔT = KomaMRICore.durs(seq)
	T0 = cumsum([0; ΔT],dims=1)
	ts = t .+ Δt
	#Eddy currents per lambda
	k = zeros(length(t), 3)
	for (i, l) in enumerate(λ)
		aux, _ =  KomaMRICore.get_eddy_currents(seq; Δt=dt, λ=l)
		k .+= α[i] .* aux
	end
	#Plot eddy currents
	p = [scatter() for j=1:4]
	p[1] = scatter(x=ts*1e3, y=k[:,1], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)", name="M2x")
	p[2] = scatter(x=ts*1e3, y=k[:,2], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)", name="M2y")
	p[3] = scatter(x=ts*1e3, y=k[:,3], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m⋅ms³)", name="M2z")
	#Layout and config
	l, config = generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0)
	plot(p, l; config)
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
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_kspace(seq)
```
"""
function plot_kspace(seq; width=nothing, height=nothing, darkmode=false)
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	#Calculations of theoretical k-space
	kspace, kspace_adc = get_kspace(seq; Δt=1) #simParams["Δt"])
	t_adc = KomaMRICore.get_adc_sampling_times(seq)
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
- `view_2d`: (`::Bool`, `=false`) boolean to use a 2D scatter plot
- `colorbar`: (`::Bool`, `=true`) boolean to show the colorbar

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
function plot_phantom_map(ph::Phantom, key::Symbol; t0=0, height=600, width=nothing, darkmode=false, view_2d=false, colorbar=true)
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
	x0 = -maximum(abs.([ph.x ph.y ph.z]))*1e2
    xf =  maximum(abs.([ph.x ph.y ph.z]))*1e2
	#Layout
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = PlotlyJS.Layout(;title=ph.name*": "*string(key),
		xaxis_title="x",
		yaxis_title="y",
		plot_bgcolor=plot_bgcolor,
		paper_bgcolor=bgcolor,
		xaxis_gridcolor=grid_color,
		yaxis_gridcolor=grid_color,
		xaxis_zerolinecolor=grid_color,
		yaxis_zerolinecolor=grid_color,
		font_color=text_color,
		scene=attr(
			xaxis=attr(title="x",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			yaxis=attr(title="y",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			zaxis=attr(title="z",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			aspectmode="manual",
			aspectratio=attr(x=1,y=1,z=1)),
		margin=attr(t=50,l=0,r=0),
		modebar=attr(orientation="h",bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
		xaxis=attr(constrain="domain"),
		yaxis=attr(scaleanchor="x"),
		hovermode="closest")
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	if view_2d
	h = PlotlyJS.scatter( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
						y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
						mode="markers",
						marker=attr(color=getproperty(ph,key)*factor,
									showscale=colorbar,
									colorscale=colormap,
									colorbar=attr(ticksuffix=unit, title=string(key)),
									cmin=cmin_key,
									cmax=cmax_key,
									size=4
									),
						text=round.(getproperty(ph,key)*factor,digits=4),
						hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
	else
	h = PlotlyJS.scatter3d( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							z=(ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							mode="markers",
							marker=attr(color=getproperty(ph,key)*factor,
										showscale=colorbar,
										colorscale=colormap,
										colorbar=attr(ticksuffix=unit, title=string(key)),
										cmin=cmin_key,
										cmax=cmax_key,
										size=2
										),
							text=round.(getproperty(ph,key)*factor,digits=4),
							hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br>z: %{z:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
	end
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
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq");

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> raw = simulate(obj, seq, sys)

julia> plot_signal(raw)
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
		if N != 1
			append!(t, t0.+(0:dt:dt*(N-1)))
		else
			append!(t, t0)
		end
		append!(signal, p.data[:,1]) #Just one coil
		#To generate gap
		append!(t, t[end])
		append!(signal, [Inf + Inf*1im])
	end
	#Show simulation blocks
	shapes = []
	annotations = []
	type_names = ["precession", "excitation"]
	if !not_Koma && show_sim_blocks
		t_sim_parts = raw.params["userParameters"]["t_sim_parts"]
		type_sim_parts = raw.params["userParameters"]["type_sim_parts"]

		current_type = -1
		for i = eachindex(t_sim_parts[1:end-1])
			aux = rect(
				xref="x", yref="paper",
				x0=t_sim_parts[i]*1e3, y0=0,
				x1=t_sim_parts[i+1]*1e3, y1=1,
				fillcolor=type_sim_parts[i] ? "Purple" : "Blue",
				opacity=.1,
				layer="below", line_width=2,
				)
			push!(shapes, aux)

			if type_sim_parts[i] != current_type
				aux = attr(xref="x", yref="paper", x=t_sim_parts[i]*1e3, y=1, showarrow=false, text=type_names[type_sim_parts[i]+1])
				push!(annotations, aux)
				current_type = type_sim_parts[i]
			end
		end
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
			annotations = annotations,
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
	PlotlyJS.plot(p, l; config)
end

"""
    str = plot_dict(dict::Dict)

Generates an HTML table based on the dictionary `dict`.

# Arguments
- `dict`: (`::Dict`) dictionary

# Returns
- `str`: (`::String`) dictionary as an HTML table
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

function plot_seqd(seq::Sequence; simParams=KomaMRICore.default_sim_params())
	seqd = KomaMRICore.discretize(seq; simParams)
	Gx = scatter(x=seqd.t*1e3, y=seqd.Gx*1e3, name="Gx", mode="markers+lines", marker_symbol=:circle)
	Gy = scatter(x=seqd.t*1e3, y=seqd.Gy*1e3, name="Gy", mode="markers+lines", marker_symbol=:circle)
	Gz = scatter(x=seqd.t*1e3, y=seqd.Gz*1e3, name="Gz", mode="markers+lines", marker_symbol=:circle)
	B1 = scatter(x=seqd.t*1e3, y=abs.(seqd.B1*1e6), name="|B1|", mode="markers+lines", marker_symbol=:circle)
	ADC = scatter(x=seqd.t[seqd.ADC]*1e3, y=zeros(sum(seqd.ADC)), name="ADC", mode="markers", marker_symbol=:x)
	plot([Gx,Gy,Gz,B1,ADC]), seqd
end
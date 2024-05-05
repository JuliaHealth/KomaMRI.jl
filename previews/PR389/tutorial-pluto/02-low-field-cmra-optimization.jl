### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> image = "https://upload.wikimedia.org/wikipedia/commons/a/a0/Textformatting.svg"
#> title = "Low-Field CMRA Optimization"
#> tags = ["CMRA", "Low Field", "Optimization"]
#> date = "2024-04-16"
#> description = "Optimizing sequence to improve SNR and fat supression."
#>
#>     [[frontmatter.author]]
#>     name = "Carlos Castillo Passi"
#>     url = "https://avatars.githubusercontent.com/u/5957134?s=400&u=fe62a2a899ced18e8b882cebde6b1eefe6a1222c&v=4"

using Markdown
using InteractiveUtils

# ╔═╡ 7deadd58-b202-4508-b4c7-686f742cb713
begin
	begin
		using Pkg
		begin
		  println("OS $(Base.Sys.MACHINE)")    # OS
		  println("Julia $VERSION")            # Julia version
		  # Koma sub-packages
		  for (_, pkg) in filter(((_, pkg),) -> occursin("KomaMRI", pkg.name), Pkg.dependencies())
		    println("$(pkg.name) $(pkg.version)")
		  end
		end
	end
end

# ╔═╡ 0b7a405e-bbb5-11ee-05ca-4b1c8567398d
using KomaMRI, PlutoPlotly # Essentials

# ╔═╡ 70dbc2bd-8b93-471d-8340-04d98a008ca6
using Suppressor, PlutoUI, ProgressLogging # Extras

# ╔═╡ a99c0c47-1b70-4362-a2f6-a7e3259606fa
md"""# Low-Field CMRA Optimization

This notebook reproduces the simulation experiments performed in our manuscript:

> **"Highly Efficient iNAV-based 3D Whole-Heart CMRA at 0.55T"**
>
> **Carlos Castillo-Passi**, Karl P. Kunze, Michael G. Crabb, Camila Muñoz, Anastasia Fotaki, Radhouene Neji, Pablo Irarrazaval, Claudia Prieto, and René M. Botnar
>
> (2024)

Submitted to Magnetic Resonance in Medicine (MRM).
"""

# ╔═╡ 1bb3e49b-1a19-4343-ac09-fbaf1cae4ba3
TableOfContents()

# ╔═╡ fc58640b-c44c-4c25-a4fd-5a0e17d7becd
md"# 1. Simulation setup"

# ╔═╡ dbf16676-64f2-4d9b-bf1b-8de06b048602
md"## 1.1. Loading required packages"

# ╔═╡ f49655cc-460e-4981-92ea-dfd6147308bf
md"Bloch simulations were performed using **KomaMRI.jl** to optimize the proposed whole-heart CMRA parameters."

# ╔═╡ 77153e5c-71bd-42e3-bae9-e4811ffa7a3d
md"""## 1.2. Scanner

We start by defining the hardware characteristics. The `sys.B0` will be used to calculate the off-resonance of fat."""

# ╔═╡ 9e397426-b60b-4b98-be8b-f7f128621c44
begin
	sys = Scanner()
	sys.B0 = 0.55
	sys.Gmax = 40.0e-3
	sys.Smax = 25.0
	sys
end

# ╔═╡ 92194fcb-582a-49ce-aad7-20b0145d40d3
md"## 1.3. Sequence"

# ╔═╡ ae72ffc5-7f8e-4907-a99e-8ad7cb8fddab
md"""
The CMRA sequence (`CMRA_iNAV_bSSFP_cardiac`) consists of:
"""

# ╔═╡ cb659118-9f22-43c3-801d-49241dee4df6
begin
	# General sequence parameters
	Trf = 500e-6  			# 500 [ms]
	B1 = 1 / (360*γ*Trf)    # B1 amplitude [uT]
	Tadc = 1e-6 			# 1us

	# Prepulses
	Tfatsat = 26.624e-3     # 26.6 [ms]
	T2prep_duration = 50e-3 # 50 [ms]

	# Acquisition
	RR = 1.0 				# 1 [s]
	dummy_heart_beats = 3 	# Steady-state
	TR = 5.3e-3             # 5.3 [ms] RF Low SAR
	TE = TR / 2 			# bSSFP condition
	iNAV_lines = 6          # FatSat-Acq delay: iNAV_lines * TR
	iNAV_flip_angle = 3.2   # 3.2 [deg]
	im_segments = 20        # Acquisitino window: im_segments * TR

	# To be optimized
	im_flip_angle = 110    # 110 [deg]
	FatSat_flip_angle = 180 # 180 [deg]

	seq_params = (;
		dummy_heart_beats,
		iNAV_lines,
		im_segments,
		iNAV_flip_angle,
		im_flip_angle,
		T2prep_duration,
		FatSat_flip_angle,
		RR
	)

	seq_params
end

# ╔═╡ a15d6b64-f8ee-4ee4-812c-d49cf5ea784d
md"""
## 1.4. Phantom
Each tissue was represented with 200 isochromats distributed along the $z$-axis to simulate gradient spoiling effects. The isochromats for each tissue were inside a 1D voxel of size $1.5\,\mathrm{mm}$. The values for $T_1$ and $T_2$ for blood, myocardial muscle, and fat at 0.55T were obtained from the work of Campbell-Washburn, et al. Fat spins were simulated using a chemical shift of $-3.4\,\mathrm{ppm}$, simulating regular fat with $T_1=183\,\mathrm{ms}$, and fast-recovering fat with $T_1=130\,\mathrm{ms}$.
"""

# ╔═╡ f0a81c9f-5616-4663-948f-a4084e1719af
begin
    fat_ppm = -3.4e-6 			# -3.4ppm fat-water frequency shift
    Niso = 200        			# 200 isochromats in spoiler direction
    Δx_voxel = 1.5e-3 			# 1.5 [mm]
    fat_freq = γ*sys.B0*fat_ppm # -80 [Hz]
    dx = Array(range(-Δx_voxel/2, Δx_voxel/2, Niso))
	md"- Phantom parameters (show/hide code)"
end

# ╔═╡ 6b870443-7be5-4287-b957-ca5c14eda89c
begin
	function FatSat(α, Δf; sample=false)
	    # FatSat design
		# cutoff_freq = sqrt(log(2) / 2) / a where B1(t) = exp(-(π t / a)^2)
		cutoff = fat_freq / π 			      # cutoff [Hz] => ≈1/10 RF power to water
		a = sqrt(log(2) / 2) / cutoff         # a [s]
		τ = range(-Tfatsat/2, Tfatsat/2, 64) # time [s]
		gauss_pulse = exp.(-(π * τ / a) .^ 2) # B1(t) [T]
		# FatSat prepulse
	    seq = Sequence()
	    seq += Grad(-8e-3, 3000e-6, 500e-6) #Spoiler1
	    seq += RF(gauss_pulse, Tfatsat, Δf)
	    α_ref = get_flip_angles(seq)[2]
	    seq *= (α/α_ref+0im)
	    if sample
	        seq += ADC(1, 1e-6)
	    end
	    seq += Grad(8e-3, 3000e-6, 500e-6) #Spoiler2
	    if sample
	        seq += ADC(1, 1e-6)
	    end
	    return seq
	end

	function T2prep(TE; sample=false)
	    seq = Sequence()
	    seq += RF(90 * B1, Trf)
	    seq += sample ? ADC(20, TE/2 - 1.5Trf) : Delay(TE/2 - 1.5Trf)
	    seq += RF(180im * B1 / 2, Trf*2)
	    seq += sample ? ADC(20, TE/2 - 1.5Trf) : Delay(TE/2 - 1.5Trf)
	    seq += RF(-90 * B1, Trf)
	    seq += Grad(8e-3, 6000e-6, 600e-6)
	    if sample
	        seq += ADC(1, 1e-6)
	    end
	    return seq
	end

	function bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample=false)
	    k = 0
	    seq = Sequence()
	    for i = 0 : iNAV_lines + im_segments - 1
	        if iNAV_lines != 0
	            m = (im_flip_angle - iNAV_flip_angle) / iNAV_lines
	            α = min( m * i + iNAV_flip_angle, im_flip_angle ) * (-1)^k
	        else
	            α = im_flip_angle * (-1)^k
	        end
	        seq += RF(α * B1, Trf)
	        if i < iNAV_lines && !sample
	            seq += Delay(TR - Trf)
	        else
	            seq += Delay(TE - Trf/2 - Tadc/2)
	            seq += ADC(1, Tadc)
	            seq += Delay(TR - TE - Tadc/2 - Trf/2)
	        end
	        k += 1
	    end
	    return seq
	end

	md"- Sequence building blocks: `T2prep`, `FatSat`, `bSSFP` (show/hide code)"
end

# ╔═╡ 7890f81e-cb15-48d2-a80c-9d73f9516056
begin
	function CMRA(
				dummy_heart_beats,
				iNAV_lines,
				im_segments,
				iNAV_flip_angle,
				im_flip_angle,
				T2prep_duration=50e-3,
				FatSat_flip_angle=180,
				RR=1.0;
				sample_recovery=zeros(Bool, dummy_heart_beats+1)
				)
		# Seq init
	    seq = Sequence()
	    for hb = 1 : dummy_heart_beats + 1
			sample = sample_recovery[hb] # Sampling recovery curve for hb
			# Generating seq blocks
	        t2p = T2prep(T2prep_duration; sample)
	        fatsat = FatSat(FatSat_flip_angle, fat_freq; sample)
	        bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle; sample)
	        # Concatenating seq blocks
	        seq += t2p
	        seq += fatsat
	        seq += bssfp
			# RR interval consideration
			RRdelay = RR  - dur(bssfp) - dur(t2p) - dur(fatsat)
	        seq += sample ? ADC(80, RRdelay) : Delay(RRdelay)
	    end
	    return seq
	end

	md"""- `CMRA` (show/hide code)

	```julia
	# Seq init
	seq = Sequence()
	for hb = 1 : dummy_heart_beats + 1
		# Generating seq blocks
		t2p = T2prep(T2prep_duration)
		fatsat = FatSat(FatSat_flip_angle, fat_freq)
		bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle, im_flip_angle)
		# Concatenating seq blocks
		seq += t2p
		seq += fatsat
		seq += bssfp
		# RR interval consideration
		RRdelay = RR  - dur(bssfp) - dur(t2p) - dur(fatsat)
		seq += Delay(RRdelay)
	end
	```"""
end

# ╔═╡ f57a2b6c-eb4c-45bd-8058-4a60b038925d
begin
	function cardiac_phantom(off; off_fat=fat_freq)
	    myocard = Phantom{Float64}(x=dx, ρ=0.6*ones(Niso), T1=701e-3*ones(Niso),
	                               T2=58e-3*ones(Niso),    Δw=2π*off*ones(Niso))
	    blood =   Phantom{Float64}(x=dx, ρ=0.7*ones(Niso), T1=1122e-3*ones(Niso),
	                               T2=263e-3*ones(Niso),   Δw=2π*off*ones(Niso))
	    fat1 =    Phantom{Float64}(x=dx, ρ=1.0*ones(Niso), T1=183e-3*ones(Niso),
	                               T2=93e-3*ones(Niso),    Δw=2π*(off_fat + off)*ones(Niso))
	    fat2 =    Phantom{Float64}(x=dx, ρ=1.0*ones(Niso), T1=130e-3*ones(Niso),
	                               T2=93e-3*ones(Niso),    Δw=2π*(off_fat + off)*ones(Niso))
	    obj = myocard + blood + fat1 + fat2
	    return obj
	end
	md"- Cardiac phantom (show/hide code)"
end

# ╔═╡ f21e9e59-25c3-4f06-8de4-792cb305eb01
md"""# 2. Simulation

Two simulation experiments were performed to optimize the sequence parameters, (1) to optimize the imaging flip angle, and (2) to optimize the FatSat flip angle.

"""

# ╔═╡ eceb326a-cab6-465e-8e5c-e835881bd3b0
md"""
## 2.0. Magnetization dynamics
"""

# ╔═╡ d54d6807-444f-4e0e-8fd6-84457974115a
md"Here we show the magnetization dynamics of the myocardium, blood, and fat signals at 0.55T."

# ╔═╡ d05dcba7-2f42-47bf-a172-6123d0113b3f
sim_params = Dict{String,Any}(
	"return_type"=>"mat",
	"sim_method"=>BlochDict(save_Mz=true),
	"Δt_rf"=>Trf,
	"gpu"=>false,
	"Nthreads"=>1
)

# ╔═╡ 37f7fd7f-5cb1-48b5-b877-b2bc23a1e7dd
begin
    seq = CMRA(seq_params...; sample_recovery=ones(Bool, dummy_heart_beats+1))
	obj = cardiac_phantom(0)
    magnetization = @suppress simulate(obj, seq, sys; sim_params)
	nothing # hide output
end

# ╔═╡ 0b6c1f72-b040-483c-969b-88bfe09b32c3
plot_seq(seq; show_adc=true, range=[2900, 3325], slider=true)

# ╔═╡ 1a62ae71-58db-49ea-ae6a-9aea66145963
begin
	phantom_T1 = plot(
		scatter(
			x=obj.x * 1e3,
			y=obj.T1 * 1e3,
			mode="markers",
			marker=attr(;
				color=obj.T1 * 1e3,
				colorscale=[
					[0.0, "black"],
					[183.0/maximum(obj.T1 .* 1e3), "green"],
					[701.0/maximum(obj.T1 .* 1e3), "blue"],
					[1122.0/maximum(obj.T1 .* 1e3), "red"],
				],
				cmin=0.0,
				cmax=1122.0,
				colorbar=attr(;ticksuffix="ms", title="T1"),
				showscale=false
			),
			showlegend=false
		)
	)
	relayout!(
		phantom_T1,
		yaxis_title="T1 [ms]",
		xaxis_title="x [mm]",
		xaxis_tickmode="array",
		xaxis_tickvals=[-1.5/2, 0.0, 1.5/2],
		yaxis_tickmode="array",
		yaxis_tickvals=unique(obj.T1 * 1e3),
		xaxis_range=[-1.5, 1.5],
		yaxis_range=[0.0, 1200.0],
		title="T1 map of 1D Phantom"
	)
	phantom_T2 = plot(
		scatter(
			x=obj.x * 1e3,
			y=obj.T2 * 1e3,
			mode="markers",
			marker=attr(;
				color=obj.T2 * 1e3,
				colorscale=[
					[0.0, "black"],
					[58.0/maximum(obj.T2 .* 1e3), "blue"],
					[93.0/maximum(obj.T2 .* 1e3), "green"],
					[263.0/maximum(obj.T2 .* 1e3), "red"],
				],
				cmin=0.0,
				cmax=263.0,
				colorbar=attr(;ticksuffix="ms", title="T2"),
				showscale=false
			),
			showlegend=false
		)
	)
	relayout!(
		phantom_T2,
		yaxis_title="T2 [ms]",
		xaxis_title="x [mm]",
		xaxis_tickmode="array",
		xaxis_tickvals=[-1.5/2, 0.0, 1.5/2],
		yaxis_tickmode="array",
		yaxis_tickvals=unique(obj.T2 * 1e3),
		xaxis_range=[-1.5, 1.5],
		yaxis_range=[0.0, 300.0],
		title="T2 map of 1D Phantom"
	)
	[phantom_T1 phantom_T2]
end

# ╔═╡ d9715bc1-49cd-4df8-8dbf-c06de42ad550
begin
    # Prep plots
    labs = ["Myocardium", "Blood", "Fat"]
	cols = ["blue", "red", "green"]
    spin_group = [(1:Niso)', (Niso+1:2Niso)', (2Niso+1:3Niso)']
    t = KomaMRICore.get_adc_sampling_times(seq)
    Mxy(i) = abs.(sum(magnetization[:,spin_group[i],1,1][:,1,:],dims=2)[:]/length(spin_group[i]))
    Mz(i) = real.(sum(magnetization[:,spin_group[i],2,1][:,1,:],dims=2)[:]/length(spin_group[i]))

    # Plot
    p0 = make_subplots(
		rows=2,
		cols=1,
		subplot_titles=["Mxy" "Mz" "Sequence"],
		shared_xaxes=true,
		vertical_spacing=0.1
	)
    for i=eachindex(spin_group)
        p1 = scatter(
			x=t, y=Mxy(i),
			name=labs[i],
			legendgroup=labs[i],
			marker_color=cols[i]
		)
        p2 = scatter(
			x=t,
			y=Mz(i),
			name=labs[i],
			legendgroup=labs[i],
			showlegend=false,
			marker_color=cols[i]
		)
        add_trace!(p0, p1, row=1, col=1)
        add_trace!(p0, p2, row=2, col=1)
    end
	relayout!(p0,
		yaxis_range=[0, 0.4],
		xaxis_range=[RR*dummy_heart_beats, RR*dummy_heart_beats+.250]
	)
    p0
end

# ╔═╡ 8dd704a4-bf50-4ddc-a832-d074bd52ad01
md"Three heartbeats were considered to achieve steady-state and the fourth heartbeat was used to measure the magnetization results in the next sections."

# ╔═╡ 39f44025-2974-4c4c-b0c2-e21399bbdb1f
md"""
## 2.1. Flip angle optimization

For the first simulation experiment, SNR was maximized by varying the imaging flip angle (between 20 and 180 deg). To optimize SNR independently of heart rate, multiple heart rates (between 55 and 85 bpm) were simulated and the mean and standard deviation of the obtained blood signal were calculated.
"""

# ╔═╡ 7b9ae381-d0fc-4d1d-b1ee-ad31a6445ff6
begin
    FAs = 20:10:180 		# flip angle [deg]
	RRs = 60 ./ (55:10:85)  # RR [s]
    mag1 = zeros(ComplexF64, im_segments, Niso*4, length(FAs), length(RRs))
    @progress for (m, RR) = enumerate(RRs), (n, im_flip_angle) = enumerate(FAs)
		seq_params1 = merge(seq_params, (; im_flip_angle, RR))
		sim_params1 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq1        = CMRA(seq_params1...)
		obj1        = cardiac_phantom(0)
		magaux = @suppress simulate(obj1, seq1, sys; sim_params=sim_params1)
		mag1[:, :, n, m] .= magaux[end-im_segments+1:end, :] # Last heartbeat
    end
end

# ╔═╡ d0377f9a-680d-4501-90ca-9ea3ab681db4
md"""## 2.2. FatSat flip angle optimization

For the second simulation experiment, the fat signal was minimized by varying the FatSat flip angle (between 20 and 250 deg) using six iNAV readouts (identified experimentally to result in good fat suppression). To be robust to $B_0$ inhomogeneities, multiple simulations with tissue frequency shifts (between $-1$ and $1\,\mathrm{ppm}$, twice of what was reported by Restivo et al.) were performed, and the mean and standard deviation of the obtained fat signal were calculated."""

# ╔═╡ d750cd5a-3c90-41ea-9942-5723a21da60a
begin
    FFAs = 20:20:250 						 # flip angle [deg]
	Δfs = (-1:0.2:1) .* (γ * sys.B0 * 1e-6)  # off-resonance Δf [s]
    mag2 = zeros(ComplexF64, im_segments, Niso*4, length(FFAs), length(Δfs))
    @progress for (m, Δf) = enumerate(Δfs), (n, FatSat_flip_angle) = enumerate(FFAs)
		seq_params2 = merge(seq_params, (; FatSat_flip_angle))
		sim_params2 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq2        = CMRA(seq_params2...)
		obj2        = cardiac_phantom(Δf)
		magaux = @suppress simulate(obj2, seq2, sys; sim_params=sim_params2)
		mag2[:, :, n, m] .= magaux[end-im_segments+1:end, :] # Last heartbeat
    end
end

# ╔═╡ 8dc11175-ebc8-407a-ab0b-6d543f849a72
begin
	# Labels
	labels = ["Myocardium", "Blood", "Fat (T₁=183 ms)", "Fat (T₁=130 ms)"]
	colors = ["blue", "red", "green", "purple"]
	spins = [(1:Niso)', ((Niso + 1):(2Niso))', ((2Niso + 1):(3Niso))', ((3Niso + 1):(4Niso))']
	mean(x, dim) = sum(x; dims=dim) / size(x, dim)
	std(x, dim; mu=mean(x, dim)) = sqrt.(sum(abs.(x .- mu) .^ 2; dims=dim) / (size(x, dim) - 1))
	md"Aux functions (show/hide code)"
end

# ╔═╡ f73082ff-a6d3-41f8-8796-4114fa89d2bb
begin
	# Reducing tissues's signal
	signal_myoc = reshape(
	    mean(abs.(mean(mag1[:, spins[1], :, :], 3)), 1), length(FAs), length(RRs)
	)
	signal_bloo = reshape(
	    mean(abs.(mean(mag1[:, spins[2], :, :], 3)), 1), length(FAs), length(RRs)
	)
	diff_bloo_myoc = abs.(signal_bloo .- signal_myoc)
	# Mean
	mean_myoc = mean(signal_myoc, 2)
	mean_bloo = mean(signal_bloo, 2)
	mean_diff = mean(diff_bloo_myoc,2)
	# Std
	std_myoc  = std(signal_myoc, 2)
	std_bloo  = std(signal_bloo, 2)
	std_diff = std(diff_bloo_myoc,2)
	# Plotting results
	# Mean
	s1 = scatter(;
	    x=FAs,
	    y=mean_myoc[:],
	    name=labels[1],
	    legendgroup=labels[1],
	    line=attr(; color=colors[1]),
	)
	s2 = scatter(;
	    x=FAs,
	    y=mean_bloo[:],
	    name=labels[2],
	    legendgroup=labels[2],
	    line=attr(; color=colors[2]),
	)
	s3 = scatter(;
		x=FAs,
		y=mean_diff[:],
		name="|Blood-Myoc|",
		legendgroup="|Blood-Myoc|",
		line=attr(color=colors[4])
	)
	# Std
	s4 = scatter(;
	    x=[FAs; reverse(FAs)],
	    y=[(mean_myoc .- std_myoc)[:]; reverse((mean_myoc .+ std_myoc)[:])],
	    name=labels[1],
	    legendgroup=labels[1],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,0,255,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s5 = scatter(;
	    x=[FAs; reverse(FAs)],
	    y=[(mean_bloo .- std_bloo)[:]; reverse((mean_bloo .+ std_bloo)[:])],
	    name=labels[2],
	    legendgroup=labels[2],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(255,0,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s6 = scatter(;
		x=[FAs; reverse(FAs)],
		y=[(mean_diff .- std_diff)[:]; reverse((mean_diff .+ std_diff)[:])],
		name="|Blood-Myoc|",legendgroup="|Blood-Myoc|",
		showlegend=false,
		fill="toself",
		fillcolor="rgba(255,0,255,0.2)",
		line=attr(color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	# Plots
	fig = plot([s1, s2, s3, s4, s5, s6])
	relayout!(
	    fig;
	    yaxis=attr(; title="Signal [a.u.]", tickmode="array"),
	    xaxis=attr(;
	        title="Flip angle [deg]",
	        tickmode="array",
	        tickvals=[FAs[1], 85, 110, 130, FAs[end]],
	        constrain="domain",
	    ),
	    font=attr(; family="CMU Serif", size=16, scaleanchor="x", scaleratio=1),
	    yaxis_range=[0, 0.3],
		xaxis_range=[FAs[1], FAs[end]],
	    width=600,
	    height=400,
	    hovermode="x unified",
	)
	fig
end

# ╔═╡ 44a31057-7b34-4c80-a273-6621c0773dc7
begin
	## Calculating results
	signal_myoc2 = reshape(
	    mean(abs.(mean(mag2[:, spins[1], :, :], 3)), 1), length(FFAs), length(Δfs)
	)
	signal_bloo2 = reshape(
	    mean(abs.(mean(mag2[:, spins[2], :, :], 3)), 1), length(FFAs), length(Δfs)
	)
	signal_fat2 = reshape(
	    mean(abs.(mean(mag2[:, spins[3], :, :], 3)), 1), length(FFAs), length(Δfs)
	)
	signal_fat22 = reshape(
	    mean(abs.(mean(mag2[:, spins[4], :, :], 3)), 1), length(FFAs), length(Δfs)
	)
	mean_myoc2 = mean(signal_myoc2, 2)
	mean_bloo2 = mean(signal_bloo2, 2)
	mean_fat2  = mean(signal_fat2, 2)
	mean_fat22 = mean(signal_fat22, 2)
	std_myoc2  = std(signal_myoc2, 2)
	std_bloo2  = std(signal_bloo2, 2)
	std_fat2   = std(signal_fat2, 2)
	std_fat22  = std(signal_fat22, 2)
	# Plotting results
	# Mean
	s12 = scatter(;
	    x=FFAs,
	    y=mean_myoc2[:],
	    name=labels[1],
	    legendgroup=labels[1],
	    line=attr(; color=colors[1]),
	)
	s22 = scatter(;
	    x=FFAs,
	    y=mean_bloo2[:],
	    name=labels[2],
	    legendgroup=labels[2],
	    line=attr(; color=colors[2]),
	)
	s32 = scatter(;
	    x=FFAs,
	    y=mean_fat2[:],
	    name=labels[3],
	    legendgroup=labels[3],
	    line=attr(; color=colors[3]),
	)
	s52 = scatter(;
	    x=FFAs,
	    y=mean_fat22[:],
	    name=labels[4],
	    legendgroup=labels[4],
	    line=attr(; color=colors[3], dash="dash"),
	)
	# Std
	s42 = scatter(;
	    x=[FFAs; reverse(FFAs)],
	    y=[(mean_myoc2 .- std_myoc2)[:]; reverse((mean_myoc2 .+ std_myoc2)[:])],
	    name=labels[1],
	    legendgroup=labels[1],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,0,255,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	s62 = scatter(;
	    x=[FFAs; reverse(FFAs)],
	    y=[(mean_bloo2 .- std_bloo2)[:]; reverse((mean_bloo2 .+ std_bloo2)[:])],
	    name=labels[2],
	    legendgroup=labels[2],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(255,0,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	s72 = scatter(;
	    x=[FFAs; reverse(FFAs)],
	    y=[(mean_fat2 .- std_fat2)[:]; reverse((mean_fat2 .+ std_fat2)[:])],
	    name=labels[3],
	    legendgroup=labels[3],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,255,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	s82 = scatter(;
	    x=[FFAs; reverse(FFAs)],
	    y=[(mean_fat22 .- std_fat22)[:]; reverse((mean_fat22 .+ std_fat22)[:])],
	    name=labels[4],
	    legendgroup=labels[4],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,255,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	# Plots
	fig2 = plot([s12, s22, s32, s42, s52, s62, s72, s82])
	relayout!(
	    fig2;
	    yaxis=attr(; title="Signal [a.u.]", tickmode="array"),
	    xaxis=attr(;
	        title="FatSat flip angle [deg]",
	        tickmode="array",
	        tickvals=[FFAs[1], 130, 150, 180, FFAs[end]],
	        constrain="domain",
	    ),
	    font=attr(; family="CMU Serif", size=16, scaleanchor="x", scaleratio=1),
	    yaxis_range=[0, 0.4],
		xaxis_range=[FFAs[1], FFAs[end]],
	    width=600,
	    height=400,
	    hovermode="x unified",
	)
	fig2
end

# ╔═╡ 3d7e7d20-a77a-48b3-ad2e-6b621227be16
md"""# References
 - **Castillo-Passi C**, Coronado R, Varela-Mattatall G, Alberola-López C, Botnar R, Irarrazaval P. KomaMRI.jl: An open-source framework for general MRI simulations with GPU acceleration. Magnetic Resonance in Medicine. 2023;90(1):329-342. [doi:10.1002/mrm.29635](doi:10.1002/mrm.29635)
 - **Campbell-Washburn AE**, Ramasawmy R, Restivo MC, et al. Opportunities in Interventional and Diagnostic Imaging by Using High-Performance Low-Field-Strength MRI. Radiology. 2019;293(2):384-393. [doi:10.1148/radiol.2019190452](doi:10.1148/radiol.2019190452)
- **Restivo MC**, Ramasawmy R, Bandettini WP, Herzka DA, Campbell-Washburn AE. Efficient spiral in-out and EPI balanced steady-state free precession cine imaging using a high-performance 0.55T MRI. Magnetic Resonance in Medicine. 2020;84(5):2364-2375. [doi:10.1002/mrm.28278](doi:10.1002/mrm.28278)
"""

# ╔═╡ abea2c43-d83e-4438-8cd3-4be06b8174b3
md"""# Reproducibility

This [Pluto notebook](https://plutojl.org/) is reproducible by default, as it has an embedded `Project.toml` and `Manifest.toml`, that store the exact package versions used to create the notebook."""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
KomaMRI = "6a340f8b-2cdf-4c04-99be-4953d9b66d0a"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Suppressor = "fd094767-a336-5f1f-9728-57cf17d0bbfb"

[compat]
KomaMRI = "~0.8.1"
PlutoPlotly = "~0.4.6"
PlutoUI = "~0.7.58"
ProgressLogging = "~0.1.4"
Suppressor = "~0.2.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "b6155d453cfb88ab9deb0f95537595210af3feea"

[[deps.AMD]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "45a1272e3f809d36431e57ab22703c6896b8908f"
uuid = "14f7f29c-3bd6-536c-9a0b-7339e30b5a3e"
version = "0.5.3"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractNFFTs]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "292e21e99dedb8621c15f185b8fdb4260bb3c429"
uuid = "7f219486-4aa7-41d6-80a7-e08ef20ceed7"
version = "0.8.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "297b6b41b66ac7cbbebb4a740844310db9fd7b8c"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.1"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "c0d491ef0b135fd7d63cbc6404286bc633329425"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.36"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "44691067188f6bd1b2289552a23e4b7572f4528d"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.9.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BFloat16s]]
deps = ["LinearAlgebra", "Printf", "Random", "Test"]
git-tree-sha1 = "2c7cc21e8678eff479978a0a2ef5ce2f51b63dff"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.5.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "7aa7ad1682f3d5754e3491bb59b8103cae28e3a3"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.40"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "3e93fcd95fe8db4704e98dbda14453a0bfc6f6c3"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.3"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BasicInterpolators]]
deps = ["LinearAlgebra", "Memoize", "Random"]
git-tree-sha1 = "3f7be532673fc4a22825e7884e9e0e876236b12a"
uuid = "26cce99e-4866-4b6d-ab74-862489e035e0"
version = "0.7.1"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Blink]]
deps = ["Base64", "Distributed", "HTTP", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Pkg", "Reexport", "Sockets", "WebIO"]
git-tree-sha1 = "bc93511973d1f949d45b0ea17878e6cb0ad484a1"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.9"

[[deps.BufferedStreams]]
git-tree-sha1 = "4ae47f9a4b1dc19897d3743ff13685925c5202ec"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.2.1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.CSSUtil]]
deps = ["Colors", "JSON", "Markdown", "Measures", "WebIO"]
git-tree-sha1 = "b9fb4b464ec10e860abe251b91d4d049934f7399"
uuid = "70588ee8-6100-5070-97c1-3cb50ed05fe8"
version = "0.1.1"

[[deps.CUDA]]
deps = ["AbstractFFTs", "Adapt", "BFloat16s", "CEnum", "CUDA_Driver_jll", "CUDA_Runtime_Discovery", "CUDA_Runtime_jll", "Crayons", "DataFrames", "ExprTools", "GPUArrays", "GPUCompiler", "KernelAbstractions", "LLVM", "LLVMLoopInfo", "LazyArtifacts", "Libdl", "LinearAlgebra", "Logging", "NVTX", "Preferences", "PrettyTables", "Printf", "Random", "Random123", "RandomNumbers", "Reexport", "Requires", "SparseArrays", "StaticArrays", "Statistics"]
git-tree-sha1 = "3dcab8a2c18ca319ea15a41d90e9528b8e93894a"
uuid = "052768ef-5323-5732-b1bb-66c8b64840ba"
version = "5.3.0"
weakdeps = ["ChainRulesCore", "SpecialFunctions"]

    [deps.CUDA.extensions]
    ChainRulesCoreExt = "ChainRulesCore"
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.CUDA_Driver_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "dc172b558adbf17952001e15cf0d6364e6d78c2f"
uuid = "4ee394cb-3365-5eb0-8335-949819d2adfc"
version = "0.8.1+0"

[[deps.CUDA_Runtime_Discovery]]
deps = ["Libdl"]
git-tree-sha1 = "38f830504358e9972d2a0c3e5d51cb865e0733df"
uuid = "1af6417a-86b4-443c-805f-a4643ffb695f"
version = "0.2.4"

[[deps.CUDA_Runtime_jll]]
deps = ["Artifacts", "CUDA_Driver_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "4ca7d6d92075906c2ce871ea8bba971fff20d00c"
uuid = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"
version = "0.12.1+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "97d79461925cdb635ee32116978fc735b9463a39"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.19"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "bfe82a708416cf00b73a3198db0859c82f741558"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.10.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d3e63d9fa13f8eaa2f06f64949e2afc593ff52c2"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.10"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArrays]]
deps = ["Adapt", "GPUArraysCore", "LLVM", "LinearAlgebra", "Printf", "Random", "Reexport", "Serialization", "Statistics"]
git-tree-sha1 = "68e8ff56a4a355a85d2784b94614491f8c900cde"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "10.1.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GPUCompiler]]
deps = ["ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "Scratch", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "1600477fba37c9fc067b9be21f5e8101f24a8865"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.26.4"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "e856eef26cf5bf2b0f95f8f4fc37553c72c8641c"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.2"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "82a471768b513dc39e471540fdadc84ff80ff997"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.3+3"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "8e59b47b9dc525b70550ca082ce85bcd7f5477cd"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.5"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ca0f6bf568b4bfc807e7537f081c81e35ceca114"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.10.0+0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.Interact]]
deps = ["CSSUtil", "InteractBase", "JSON", "Knockout", "Observables", "OrderedCollections", "Reexport", "WebIO", "Widgets"]
git-tree-sha1 = "c5091992248c7134af7c90554305c600d5d9012b"
uuid = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
version = "0.10.5"

[[deps.InteractBase]]
deps = ["Base64", "CSSUtil", "Colors", "Dates", "JSExpr", "JSON", "Knockout", "Observables", "OrderedCollections", "Random", "WebIO", "Widgets"]
git-tree-sha1 = "aa5daeff326db0a9126a225b58ca04ae12f57259"
uuid = "d3863d7c-f0c8-5437-a7b4-3ae773c01009"
version = "0.10.10"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "896385798a8d49a255c398bd49162062e4a4c435"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.13"
weakdeps = ["Dates"]

    [deps.InverseFunctions.extensions]
    DatesExt = "Dates"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaNVTXCallbacks_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "af433a10f3942e882d3c671aacb203e006a5808f"
uuid = "9c1d0b0a-7046-5b2e-a33f-ea22f176ac7e"
version = "0.2.1+0"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2ef87eeaa28713cb010f9fb0be288b6c1a4ecd53"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.1.0+0"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "Requires", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "ed7167240f40e62d97c1f5f7735dea6de3cc5c49"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.18"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.Knockout]]
deps = ["JSExpr", "JSON", "Observables", "Test", "WebIO"]
git-tree-sha1 = "91835de56d816864f1c38fb5e3fad6eb1e741271"
uuid = "bcebb21b-c2e3-54f8-a781-646b90f6d2cc"
version = "0.2.6"

[[deps.KomaMRI]]
deps = ["AssetRegistry", "Blink", "FFTW", "Interact", "KomaMRICore", "KomaMRIFiles", "KomaMRIPlots", "MAT", "MRIReco", "Pkg", "Reexport"]
git-tree-sha1 = "5580dcf41811d073c0e69109eabec8e9fcf4e909"
uuid = "6a340f8b-2cdf-4c04-99be-4953d9b66d0a"
version = "0.8.1"

[[deps.KomaMRIBase]]
deps = ["Interpolations", "MAT", "MRIBase", "Parameters", "Pkg", "Reexport"]
git-tree-sha1 = "7d2cb6ae57f90d06483121ca2d3f1d161dd25fe4"
uuid = "d0bc0b20-b151-4d03-b2a4-6ca51751cb9c"
version = "0.8.3"

[[deps.KomaMRICore]]
deps = ["Adapt", "CUDA", "Functors", "KomaMRIBase", "Pkg", "ProgressMeter", "Reexport", "ThreadsX"]
git-tree-sha1 = "32b3b97e7e9c9250647b0304403592646159e308"
uuid = "4baa4f4d-2ae9-40db-8331-a7d1080e3f4e"
version = "0.8.2"

[[deps.KomaMRIFiles]]
deps = ["FileIO", "HDF5", "KomaMRIBase", "MAT", "MRIFiles", "Pkg", "Reexport", "Scanf"]
git-tree-sha1 = "9720fd5fd68d0ee5d805a60cc7ed8e6306dc8988"
uuid = "fcf631a6-1c7e-4e88-9e64-b8888386d9dc"
version = "0.8.2"

[[deps.KomaMRIPlots]]
deps = ["Interpolations", "Kaleido_jll", "KomaMRIBase", "MAT", "Pkg", "PlotlyJS", "Reexport"]
git-tree-sha1 = "103092ea5e15c09ed65b6f9d891bf724e53c262b"
uuid = "76db0263-63f3-4d26-bb9a-5dba378db904"
version = "0.8.2"
weakdeps = ["PlutoPlotly"]

    [deps.KomaMRIPlots.extensions]
    KomaPlotsPlutoPlotlyExt = "PlutoPlotly"

[[deps.LDLFactorizations]]
deps = ["AMD", "LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "70f582b446a1c3ad82cf87e62b878668beef9d13"
uuid = "40e66cde-538c-5869-a4ad-c39174c6795b"
version = "0.10.1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Requires", "Unicode"]
git-tree-sha1 = "839c82932db86740ae729779e610f07a1640be9a"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "6.6.3"
weakdeps = ["BFloat16s"]

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "88b916503aac4fb7f701bb625cd84ca5dd1677bc"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.29+0"

[[deps.LLVMLoopInfo]]
git-tree-sha1 = "2e5c102cfc41f48ae4740c7eca7743cc7e7b75ea"
uuid = "8b046642-f1f6-4319-8d3c-209ddc03c586"
version = "1.0.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearOperatorCollection]]
deps = ["InteractiveUtils", "LinearAlgebra", "LinearOperators", "Random", "Reexport", "SparseArrays"]
git-tree-sha1 = "9bec337396ad9c46fa9f662606718755dbcda5eb"
uuid = "a4a2c56f-fead-462a-a3ab-85921a5f2575"
version = "1.2.1"
weakdeps = ["FFTW", "NFFT", "Wavelets"]

    [deps.LinearOperatorCollection.extensions]
    LinearOperatorFFTWExt = "FFTW"
    LinearOperatorNFFTExt = ["NFFT", "FFTW"]
    LinearOperatorWaveletExt = "Wavelets"

[[deps.LinearOperators]]
deps = ["FastClosures", "LDLFactorizations", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "TimerOutputs"]
git-tree-sha1 = "f06df3a46255879cbccae1b5b6dcb16994c31be7"
uuid = "5c8ed15e-5a4c-59e4-a42b-c7e8811fb125"
version = "2.7.0"
weakdeps = ["ChainRulesCore"]

    [deps.LinearOperators.extensions]
    LinearOperatorsChainRulesCoreExt = "ChainRulesCore"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LowRankApprox]]
deps = ["FFTW", "LinearAlgebra", "LowRankMatrices", "Nullables", "Random", "SparseArrays"]
git-tree-sha1 = "031af63ba945e23424815014ba0e59c28f5aed32"
uuid = "898213cb-b102-5a47-900c-97e73b919f73"
version = "0.5.5"

[[deps.LowRankMatrices]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7c8664b2f3d5c3d9b77605c03d53b18813e79b0f"
uuid = "e65ccdef-c354-471a-8090-89bec1c20ec3"
version = "1.0.1"
weakdeps = ["FillArrays"]

    [deps.LowRankMatrices.extensions]
    LowRankMatricesFillArraysExt = "FillArrays"

[[deps.MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "ed1cf0a322d78cee07718bed5fd945e2218c35a1"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.6"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "656036b9ed6f942d35e536e249600bc31d0f9df8"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.2.0+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8f6af051b9e8ec597fa09d8885ed79fd582f33c9"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.10"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "77c3bd69fdb024d75af38713e883d0f249ce19c2"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.3.2+0"

[[deps.MRIBase]]
deps = ["AbstractNFFTs", "LinearAlgebra", "NFFTTools"]
git-tree-sha1 = "a96ba4ecb38712ebc395891e375dac99adb8bd72"
uuid = "f7771a9a-6e57-4e71-863b-6e4b6a2f17df"
version = "0.4.2"

[[deps.MRIFiles]]
deps = ["FileIO", "HDF5", "LightXML", "LinearAlgebra", "MRIBase", "Printf", "Reexport"]
git-tree-sha1 = "afe2461175ccfcaff74ac211f9290e990adb6262"
uuid = "5a6f062f-bf45-497d-b654-ad17aae2a530"
version = "0.3.0"

[[deps.MRIOperators]]
deps = ["Distributions", "FLoops", "LinearAlgebra", "LinearOperatorCollection", "LinearOperators", "LowRankApprox", "MRIBase", "NFFT", "Printf", "Reexport", "StatsBase", "Wavelets"]
git-tree-sha1 = "ad33e3745db181e93d0ba3b0309ef959ec5e517c"
uuid = "fb1137e3-90a6-46ce-a672-6e1e53d120f2"
version = "0.2.1"

[[deps.MRIReco]]
deps = ["AxisArrays", "FLoops", "LinearAlgebra", "MRIBase", "MRIOperators", "PrecompileTools", "ProgressMeter", "Random", "Reexport", "RegularizedLeastSquares", "Unitful"]
git-tree-sha1 = "648ebc967ac4d42fdd57c12a9918a1119a7cbf9a"
uuid = "bdf86e05-2d2b-5731-a332-f3fe1f9e047f"
version = "0.8.1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f12a29c4400ba812841c6ace3f4efbb6dbb3ba01"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "a7cefa21a2ff993bff0456bf7521f46fc077ddf1"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.19"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "MbedTLS", "Pkg", "Sockets"]
git-tree-sha1 = "7295d849103ac4fcbe3b2e439f229c5cc77b9b69"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "1.0.2"

[[deps.NFFT]]
deps = ["AbstractNFFTs", "BasicInterpolators", "Distributed", "FFTW", "FLoops", "LinearAlgebra", "Printf", "Random", "Reexport", "SnoopPrecompile", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "93a5f32dd6cf09456b0b81afcb8fc29f06535ffd"
uuid = "efe261a4-0d2b-5849-be55-fc731d526b0d"
version = "0.13.3"

[[deps.NFFTTools]]
deps = ["AbstractFFTs", "AbstractNFFTs", "FFTW", "LinearAlgebra"]
git-tree-sha1 = "d6a68b7ffbd50b4c99e514a1a6fb8ce84f6e247e"
uuid = "7424e34d-94f7-41d6-98a0-85abaf1b6c91"
version = "0.2.6"

[[deps.NVTX]]
deps = ["Colors", "JuliaNVTXCallbacks_jll", "Libdl", "NVTX_jll"]
git-tree-sha1 = "53046f0483375e3ed78e49190f1154fa0a4083a1"
uuid = "5da4648a-3479-48b8-97b9-01cb529c0a1f"
version = "0.3.4"

[[deps.NVTX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ce3269ed42816bf18d500c9f63418d4b0d9f5a3b"
uuid = "e98f9f5b-d649-5603-91fd-7774390e6439"
version = "3.1.0+2"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Nullables]]
git-tree-sha1 = "8f87854cc8f3685a60689d8edecaa29d2251979b"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "1.0.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "e64b4f5ea6b7389f6f046d13d4896a8f9c1ba71e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e25c1778a98e34219a00455d6e4384e017ea9762"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.6+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "PlotlyKaleido", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "e62d886d33b81c371c9d4e2f70663c0637f19459"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.13"

    [deps.PlotlyJS.extensions]
    CSVExt = "CSV"
    DataFramesExt = ["DataFrames", "CSV"]
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyJS.weakdeps]
    CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlotlyKaleido]]
deps = ["Base64", "JSON", "Kaleido_jll"]
git-tree-sha1 = "2650cd8fb83f73394996d507b3411a7316f6f184"
uuid = "f2990250-8cf9-495f-b13a-cce12b45703c"
version = "2.2.4"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "1ae939782a5ce9a004484eab5416411c7190d3ce"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.4.6"
weakdeps = ["PlotlyKaleido", "Unitful"]

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "81a2a9462003a423fdc59e2a3ff84cde93c4637b"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.7"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "763a8ceb07833dd51bb9e3bbca372de32c0605ad"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "02d31ad62838181c1a3a5fd23a1ce5914a643601"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.3"

[[deps.RegularizedLeastSquares]]
deps = ["FLoops", "FastClosures", "IterativeSolvers", "LinearAlgebra", "LinearOperators", "ProgressMeter", "Random", "SparseArrays", "StatsBase", "VectorizationBase"]
git-tree-sha1 = "d76fe5bbb6191f84ac2cc88eb67471800743c3aa"
uuid = "1e9c538a-f78c-5de5-8ffb-0b6dbe892d23"
version = "0.10.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.Scanf]]
deps = ["BufferedStreams"]
git-tree-sha1 = "cd8f9da92f751a1fc421120c66971d008d695cdb"
uuid = "6ef1bc8b-493b-44e1-8d40-549aa65c4b41"
version = "0.5.4"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d2fdac9ff3906e27f7a618d47b676941baa6c80c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.10"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Suppressor]]
deps = ["Logging"]
git-tree-sha1 = "9143c41bd539a8885c79728b9dedb0ce47dc9819"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.7"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["Accessors", "ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "SplittablesBase", "Transducers"]
git-tree-sha1 = "70bd8244f4834d46c3d68bd09e7792d8f571ef04"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.12"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
git-tree-sha1 = "71509f04d045ec714c4748c785a59045c3736349"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.7"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "3064e780dbb8a9296ebb3af8f440f787bb5332af"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.80"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "323e3d0acf5e78a56dfae7bd8928c989b4f3083e"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.3"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "ac377f0a248753a1b1d58bbc92a64f5a726dfb71"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.66"

[[deps.Wavelets]]
deps = ["DSP", "FFTW", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f514f9b16f6a15552c6aad7b03afc7b9a8478ef4"
uuid = "29a6e085-ba6d-5f35-a997-948ac2efa89a"
version = "0.10.0"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "4162e95e05e79922e44b9952ccbc262832e4ad07"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.6.0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "532e22cf7be8462035d092ff21fada7527e2c488"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.6+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "46bf7be2917b59b761247be3f317ddf75e50e997"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─a99c0c47-1b70-4362-a2f6-a7e3259606fa
# ╟─1bb3e49b-1a19-4343-ac09-fbaf1cae4ba3
# ╟─fc58640b-c44c-4c25-a4fd-5a0e17d7becd
# ╟─dbf16676-64f2-4d9b-bf1b-8de06b048602
# ╟─f49655cc-460e-4981-92ea-dfd6147308bf
# ╠═0b7a405e-bbb5-11ee-05ca-4b1c8567398d
# ╠═70dbc2bd-8b93-471d-8340-04d98a008ca6
# ╟─77153e5c-71bd-42e3-bae9-e4811ffa7a3d
# ╟─9e397426-b60b-4b98-be8b-f7f128621c44
# ╟─92194fcb-582a-49ce-aad7-20b0145d40d3
# ╟─ae72ffc5-7f8e-4907-a99e-8ad7cb8fddab
# ╟─6b870443-7be5-4287-b957-ca5c14eda89c
# ╟─7890f81e-cb15-48d2-a80c-9d73f9516056
# ╟─cb659118-9f22-43c3-801d-49241dee4df6
# ╠═0b6c1f72-b040-483c-969b-88bfe09b32c3
# ╟─a15d6b64-f8ee-4ee4-812c-d49cf5ea784d
# ╟─f0a81c9f-5616-4663-948f-a4084e1719af
# ╟─f57a2b6c-eb4c-45bd-8058-4a60b038925d
# ╟─1a62ae71-58db-49ea-ae6a-9aea66145963
# ╟─f21e9e59-25c3-4f06-8de4-792cb305eb01
# ╟─eceb326a-cab6-465e-8e5c-e835881bd3b0
# ╟─d54d6807-444f-4e0e-8fd6-84457974115a
# ╟─d05dcba7-2f42-47bf-a172-6123d0113b3f
# ╠═37f7fd7f-5cb1-48b5-b877-b2bc23a1e7dd
# ╟─d9715bc1-49cd-4df8-8dbf-c06de42ad550
# ╟─8dd704a4-bf50-4ddc-a832-d074bd52ad01
# ╟─39f44025-2974-4c4c-b0c2-e21399bbdb1f
# ╠═7b9ae381-d0fc-4d1d-b1ee-ad31a6445ff6
# ╟─f73082ff-a6d3-41f8-8796-4114fa89d2bb
# ╟─d0377f9a-680d-4501-90ca-9ea3ab681db4
# ╠═d750cd5a-3c90-41ea-9942-5723a21da60a
# ╟─44a31057-7b34-4c80-a273-6621c0773dc7
# ╟─8dc11175-ebc8-407a-ab0b-6d543f849a72
# ╟─3d7e7d20-a77a-48b3-ad2e-6b621227be16
# ╟─abea2c43-d83e-4438-8cd3-4be06b8174b3
# ╟─7deadd58-b202-4508-b4c7-686f742cb713
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

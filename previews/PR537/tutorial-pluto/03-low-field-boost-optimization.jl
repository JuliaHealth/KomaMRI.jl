### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> image = "https://upload.wikimedia.org/wikipedia/commons/a/a0/Textformatting.svg"
#> title = "Low-Field BOOST Optimization"
#> tags = ["BOOST", "Low Field", "Optimization"]
#> date = "2024-10-03"
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
using KomaMRICore, KomaMRIPlots, PlutoPlotly # Essentials

# ╔═╡ 70dbc2bd-8b93-471d-8340-04d98a008ca6
using Suppressor, PlutoUI, ProgressLogging # Extras

# ╔═╡ a99c0c47-1b70-4362-a2f6-a7e3259606fa
md"""# Low-Field BOOST Optimization

This notebook reproduces the simulation experiments performed in our manuscript:

> **"Simultaneous iNAV-Based 3D Whole-Heart Bright-Blood and Black-Blood Imaging at 0.55T"**
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
md"Bloch simulations were performed using **KomaMRI.jl** to optimize the proposed whole-heart BOOST parameters."

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
The BOOST sequence (`BOOST`) consists of:
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
	im_flip_angle = [110, 80] # 80 [deg]
	FatSat_flip_angle = 180   # 180 [deg]
	IR_inversion_time = 90e-3 # 90 [ms] 

	seq_params = (;
		dummy_heart_beats,
		iNAV_lines,
		im_segments,
		iNAV_flip_angle,
		im_flip_angle,
		T2prep_duration,
		IR_inversion_time,
		FatSat_flip_angle,
		RR
	)

	seq_params
end

# ╔═╡ 42417eba-46e0-400e-a874-299533362c41
md"Below we are showing the **bright-blood contrast** (containing a T2p-IR pulse):"

# ╔═╡ 4aab8a6f-d2ba-46f0-a1bf-ddc0dcf8437f
md"and **reference contrast** (only containing a FatSat pulse), later used to obtain the black-blood contrast by subtracting the bright-blood contrast. Note that this contrast has a different flip angle."

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
		τ = range(-Tfatsat/2, Tfatsat/2, 64)  # time [s]
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
	    seq += Grad(8e-3, 6000e-6, 600e-6) #Spoiler3
	    if sample
	        seq += ADC(1, 1e-6)
	    end
	    return seq
	end

	function IR(IR_delay; sample=false)
		# Generating HS pulse
		# Based on: https://onlinelibrary.wiley.com/doi/epdf/10.1002/jmri.26021
		# Params
		flip_angle = 900;    # Peak amplitude (deg)
		Trf = 10240e-6;      # Pulse duration (ms)
		β = 6.7e2;           # frequency modulation param (rad/s)
		μ = 5;               # phase modulation parameter (dimensionless)
		fmax = μ * β / (2π); # 2fmax = BW
		# RF pulse
		t = range(-Trf/2, Trf/2, 201);
		B1 = sech.(β .* t);
		Δf = fmax  .* tanh.(β .* t);
		# Spoiler length
		spoiler_time = 6000e-6
		spoiler_rise_fall = 600e-6
		# Prepulse
	    seq = Sequence()
	    seq += RF(B1, Trf, Δf) # FM modulated pulse
		seq = (flip_angle / get_flip_angles(seq)[1] + 0.0im) * seq # RF scaling
	    seq += Grad(8e-3, spoiler_time, spoiler_rise_fall) #Spoiler3
	    if sample
	        seq += ADC(11, IR_delay - spoiler_time - 2spoiler_rise_fall)
		else
			seq += Delay(IR_delay - spoiler_time - 2spoiler_rise_fall)
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

	md"- Sequence building blocks: `T2prep`, `IR`, `FatSat`, `bSSFP` (show/hide code)"
end

# ╔═╡ 7890f81e-cb15-48d2-a80c-9d73f9516056
begin
	function BOOST(
				dummy_heart_beats,
				iNAV_lines,
				im_segments,
				iNAV_flip_angle,
				im_flip_angle,
				T2prep_duration=50e-3,
				IR_inversion_time=90e-3,
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
			ir = IR(IR_inversion_time - iNAV_lines * TR - Trf - TE; sample)
	        fatsat = FatSat(FatSat_flip_angle, fat_freq; sample)
	        # Magnetization preparations
			for contrast = 1:2
				preps = Sequence()
				if contrast == 1 # Bright-blood contrast
		        	preps += t2p
					preps += ir
				else # Reference contrast
		        	preps += fatsat
				end
				# Contrst dependant flip angle
				bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle,
					im_flip_angle[contrast]; sample)
				# Concatenating seq blocks
				seq += preps
		        seq += bssfp
				# RR interval consideration
				RRdelay = RR  - dur(bssfp) - dur(preps)
		        seq += sample ? ADC(80, RRdelay) : Delay(RRdelay)
			end
	    end
	    return seq
	end

	md"""- `BOOST` (show/hide code)

	```julia
	# Seq init
	seq = Sequence()
	for hb = 1 : dummy_heart_beats + 1
		# Generating seq blocks
		t2p = T2prep(T2prep_duration)
		ir = IR(IR_inversion_time)
		fatsat = FatSat(FatSat_flip_angle, fat_freq)
		# Magnetization preparations
		for contrast = 1:2
			preps = Sequence()
			if contrast == 1 # Bright-blood contrast
				preps += t2p
				preps += ir
			else # Reference contrast
				preps += fatsat
			end
			# Contrst dependant flip angle
			bssfp = bSSFP(iNAV_lines, im_segments, iNAV_flip_angle,
				im_flip_angle[contrast]; sample)
			# Concatenating seq blocks
			seq += preps
			seq += bssfp
			# RR interval consideration
			RRdelay = RR  - dur(bssfp) - dur(preps)
			seq += Delay(RRdelay)
		end
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
    seq = BOOST(seq_params...; sample_recovery=ones(Bool, dummy_heart_beats+1))
	obj = cardiac_phantom(0)
    magnetization = @suppress simulate(obj, seq, sys; sim_params)
	nothing # hide output
end

# ╔═╡ 0b6c1f72-b040-483c-969b-88bfe09b32c3
plot_seq(seq; range=[5990, 6280], slider=true)

# ╔═╡ 88eb41a5-d8c2-4f0e-b379-a8b05a341a82
plot_seq(seq; range=[6900, 7190], slider=true)

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
md"Three heartbeats per contrast were considered to achieve steady-state and the fourth heartbeat for each contrast was used to measure the magnetization results in the next sections."

# ╔═╡ 39f44025-2974-4c4c-b0c2-e21399bbdb1f
md"""
## 2.1. Reference contrast SNR: Flip angle optimization

For the first simulation experiment, SNR was maximized by varying the imaging flip angle (between 20 and 180 deg). To optimize SNR independently of heart rate, multiple heart rates (between 55 and 85 bpm) were simulated and the mean and standard deviation of the obtained blood signal were calculated.
"""

# ╔═╡ 7b9ae381-d0fc-4d1d-b1ee-ad31a6445ff6
begin
    FAs = 20:10:180          # flip angle [deg]
	RRs = 60 ./ (55:10:85)  # RR [s]
    mag1 = zeros(ComplexF64, im_segments, Niso*4, length(FAs), length(RRs))
    @progress for (m, RR) = enumerate(RRs), (n, α) = enumerate(FAs)
		seq_params1 = merge(seq_params, (; im_flip_angle=[110, α], RR))
		sim_params1 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq1        = BOOST(seq_params1...)
		obj1        = cardiac_phantom(0)
		magaux = @suppress simulate(obj1, seq1, sys; sim_params=sim_params1)
		mag1[:, :, n, m] .= magaux[end-im_segments+1:end, :] # Last heartbeat
    end
end

# ╔═╡ d0377f9a-680d-4501-90ca-9ea3ab681db4
md"""## 2.2. Reference contrast fat: FatSat flip angle optimization

For the second simulation experiment, the fat signal was minimized by varying the FatSat flip angle (between 20 and 250 deg) using six iNAV readouts (identified experimentally to result in good fat suppression). To be robust to $B_0$ inhomogeneities, multiple simulations with tissue frequency shifts (between $-1$ and $1\,\mathrm{ppm}$, twice of what was reported by Restivo et al.) were performed, and the mean and standard deviation of the obtained fat signal were calculated."""

# ╔═╡ d750cd5a-3c90-41ea-9942-5723a21da60a
begin
    FFAs = 20:20:250 						 # flip angle [deg]
	Δfs = (-1:0.2:1) .* (γ * sys.B0 * 1e-6)  # off-resonance Δf [s]
    mag2 = zeros(ComplexF64, im_segments, Niso*4, length(FFAs), length(Δfs))
    @progress for (m, Δf) = enumerate(Δfs), (n, FatSat_flip_angle) = enumerate(FFAs)
		seq_params2 = merge(seq_params, (; FatSat_flip_angle))
		sim_params2 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq2        = BOOST(seq_params2...)
		obj2        = cardiac_phantom(Δf)
		magaux = @suppress simulate(obj2, seq2, sys; sim_params=sim_params2)
		mag2[:, :, n, m] .= magaux[end-im_segments+1:end, :] # Last heartbeat
    end
end

# ╔═╡ b0b53632-e6b8-47e9-8d3e-79f1e599315a
md" ## 2.3. Bright-blood SNR:"

# ╔═╡ 14cf6859-a4b9-4671-b022-659781c55144
begin
    mag4 = zeros(ComplexF64, im_segments, Niso*4, length(FAs), length(RRs))
    @progress for (m, RR) = enumerate(RRs), (n, α) = enumerate(FAs)
		seq_params4 = merge(seq_params, (; im_flip_angle=[α, 80], RR))
		sim_params4 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq4        = BOOST(seq_params4...)
		obj4        = cardiac_phantom(0)
		magaux = @suppress simulate(obj4, seq4, sys; sim_params=sim_params4)
		mag4[:, :, n, m] .= magaux[end-2im_segments+1:end-im_segments, :] # Bright-Blood
    end
end

# ╔═╡ c64d38cb-2433-4b90-abbf-1fc938f70584
md"## 2.4. Bright-Blood fat: Inversion delay"

# ╔═╡ c952ecf1-25ef-4b48-9c8f-e53ded302629
begin
    TIs = (50:10:130) # Inversion delay [ms]
    mag3 = zeros(ComplexF64, im_segments, Niso*4, length(TIs), length(Δfs))
    @progress for (m, Δf) = enumerate(Δfs), (n, TI) = enumerate(TIs)
		seq_params3 = merge(seq_params, (; IR_inversion_time=TI * 1e-3, RR))
		sim_params3 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq3        = BOOST(seq_params3...)
		obj3        = cardiac_phantom(Δf)
		magaux = @suppress simulate(obj3, seq3, sys; sim_params=sim_params3)
		mag3[:, :, n, m] .= magaux[end-2im_segments+1:end-im_segments, :] # Bright-Blood
    end
end

# ╔═╡ 85834365-238b-4193-a0c0-d1859416ab6c
md" ## 2.5. Bright-blood contrast: T2prep duration"

# ╔═╡ 50b551d8-a212-4319-b25b-1f8ecdd2a491
begin
    T2ps = (20:10:80) # T2prep duration [ms]
    mag5bb = zeros(ComplexF64, im_segments, Niso*4, length(T2ps), length(RRs))
	mag5rf = zeros(ComplexF64, im_segments, Niso*4, length(T2ps), length(RRs))
    @progress for (m, RR) = enumerate(RRs), (n, T2p) = enumerate(T2ps)
		seq_params5 = merge(seq_params, (; T2prep_duration=T2p * 1e-3, RR))
		sim_params5 = merge(sim_params, Dict("sim_method"=>BlochDict()))
		seq5        = BOOST(seq_params5...)
		obj5        = cardiac_phantom(0)
		magaux = @suppress simulate(obj5, seq5, sys; sim_params=sim_params5)
		# Bright-Blood
		mag5bb[:, :, n, m] .= magaux[end-2im_segments+1:end-im_segments, :]
		# Reference contrast
		mag5rf[:, :, n, m] .= magaux[end-im_segments+1:end, :] 
    end
end

# ╔═╡ 4b1dc8bf-feb1-42b4-94f7-1b27975d944c
md"## 2.6. Black-blood contrast: T2prep duration"

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
	        tickvals=[FAs[1], 80, 110, 130, FAs[end]],
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

# ╔═╡ 6649f46e-3565-4cd6-86a5-9b03d00cc3db
begin
	# Reducing tissues's signal
	signal_myoc4 = reshape(
	    mean(abs.(mean(mag4[:, spins[1], :, :], 3)), 1), length(FAs), length(RRs)
	)
	signal_bloo4 = reshape(
	    mean(abs.(mean(mag4[:, spins[2], :, :], 3)), 1), length(FAs), length(RRs)
	)
	diff_bloo_myoc4 = abs.(signal_bloo4 .- signal_myoc4)
	# Mean
	mean_myoc4 = mean(signal_myoc4, 2)
	mean_bloo4 = mean(signal_bloo4, 2)
	mean_diff4 = mean(diff_bloo_myoc4,2)
	# Std
	std_myoc4  = std(signal_myoc4, 2)
	std_bloo4  = std(signal_bloo4, 2)
	std_diff4 = std(diff_bloo_myoc4,2)
	# Plotting results
	# Mean
	s14 = scatter(;
	    x=FAs,
	    y=mean_myoc4[:],
	    name=labels[1],
	    legendgroup=labels[1],
	    line=attr(; color=colors[1]),
	)
	s24 = scatter(;
	    x=FAs,
	    y=mean_bloo4[:],
	    name=labels[2],
	    legendgroup=labels[2],
	    line=attr(; color=colors[2]),
	)
	s34 = scatter(;
		x=FAs,
		y=mean_diff4[:],
		name="|Blood-Myoc|",
		legendgroup="|Blood-Myoc|",
		line=attr(color=colors[4])
	)
	# Std
	s44 = scatter(;
	    x=[FAs; reverse(FAs)],
	    y=[(mean_myoc4 .- std_myoc4)[:]; reverse((mean_myoc4 .+ std_myoc4)[:])],
	    name=labels[1],
	    legendgroup=labels[1],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,0,255,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s54 = scatter(;
	    x=[FAs; reverse(FAs)],
	    y=[(mean_bloo4 .- std_bloo4)[:]; reverse((mean_bloo4 .+ std_bloo4)[:])],
	    name=labels[2],
	    legendgroup=labels[2],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(255,0,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s64 = scatter(;
		x=[FAs; reverse(FAs)],
		y=[(mean_diff4 .- std_diff4)[:]; reverse((mean_diff4 .+ std_diff4)[:])],
		name="|Blood-Myoc|",legendgroup="|Blood-Myoc|",
		showlegend=false,
		fill="toself",
		fillcolor="rgba(255,0,255,0.2)",
		line=attr(color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	# Plots
	fig4 = plot([s14, s24, s34, s44, s54, s64])
	relayout!(
	    fig4;
	    yaxis=attr(; title="Signal [a.u.]", tickmode="array"),
	    xaxis=attr(;
	        title="Flip angle [deg]",
	        tickmode="array",
	        tickvals=[FAs[1], 110, FAs[end]],
	        constrain="domain",
	    ),
	    font=attr(; family="CMU Serif", size=16, scaleanchor="x", scaleratio=1),
	    yaxis_range=[0, 0.3],
		xaxis_range=[FAs[1], FAs[end]],
	    width=600,
	    height=400,
	    hovermode="x unified",
	)
	fig4
end

# ╔═╡ a8751a97-c131-407b-b211-573660452142
begin
	## Calculating results
	signal_myoc3 = reshape(
	    mean(abs.(mean(mag3[:, spins[1], :, :], 3)), 1), length(TIs), length(Δfs)
	)
	signal_bloo3 = reshape(
	    mean(abs.(mean(mag3[:, spins[2], :, :], 3)), 1), length(TIs), length(Δfs)
	)
	signal_fat3 = reshape(
	    mean(abs.(mean(mag3[:, spins[3], :, :], 3)), 1), length(TIs), length(Δfs)
	)
	signal_fat33 = reshape(
	    mean(abs.(mean(mag3[:, spins[4], :, :], 3)), 1), length(TIs), length(Δfs)
	)
	mean_myoc3 = mean(signal_myoc3, 2)
	mean_bloo3 = mean(signal_bloo3, 2)
	mean_fat3  = mean(signal_fat3, 2)
	mean_fat33 = mean(signal_fat33, 2)
	std_myoc3  = std(signal_myoc3, 2)
	std_bloo3  = std(signal_bloo3, 2)
	std_fat3   = std(signal_fat3, 2)
	std_fat33  = std(signal_fat33, 2)
	# Plotting results
	# Mean
	s13 = scatter(;
	    x=TIs,
	    y=mean_myoc3[:],
	    name=labels[1],
	    legendgroup=labels[1],
	    line=attr(; color=colors[1]),
	)
	s23 = scatter(;
	    x=TIs,
	    y=mean_bloo3[:],
	    name=labels[2],
	    legendgroup=labels[2],
	    line=attr(; color=colors[2]),
	)
	s33 = scatter(;
	    x=TIs,
	    y=mean_fat3[:],
	    name=labels[3],
	    legendgroup=labels[3],
	    line=attr(; color=colors[3]),
	)
	s53 = scatter(;
	    x=TIs,
	    y=mean_fat33[:],
	    name=labels[4],
	    legendgroup=labels[4],
	    line=attr(; color=colors[3], dash="dash"),
	)
	# Std
	s43 = scatter(;
	    x=[TIs; reverse(TIs)],
	    y=[(mean_myoc3 .- std_myoc3)[:]; reverse((mean_myoc3 .+ std_myoc3)[:])],
	    name=labels[1],
	    legendgroup=labels[1],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,0,255,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	s63 = scatter(;
	    x=[TIs; reverse(TIs)],
	    y=[(mean_bloo3 .- std_bloo3)[:]; reverse((mean_bloo3 .+ std_bloo3)[:])],
	    name=labels[2],
	    legendgroup=labels[2],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(255,0,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	s73 = scatter(;
	    x=[TIs; reverse(TIs)],
	    y=[(mean_fat3 .- std_fat3)[:]; reverse((mean_fat3 .+ std_fat3)[:])],
	    name=labels[3],
	    legendgroup=labels[3],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,255,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	s83 = scatter(;
	    x=[TIs; reverse(TIs)],
	    y=[(mean_fat33 .- std_fat33)[:]; reverse((mean_fat33 .+ std_fat33)[:])],
	    name=labels[4],
	    legendgroup=labels[4],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,255,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none",
	)
	# Plots
	fig3 = plot([s13, s23, s33, s43, s53, s63, s73, s83])
	relayout!(
	    fig3;
	    yaxis=attr(; title="Signal [a.u.]", tickmode="array"),
	    xaxis=attr(;
	        title="Inversion Time [ms]",
	        tickmode="array",
	        tickvals=[TIs[1], 70, 90, TIs[end]],
	        constrain="domain",
	    ),
	    font=attr(; family="CMU Serif", size=16, scaleanchor="x", scaleratio=1),
	    yaxis_range=[0, 0.3],
		xaxis_range=[TIs[1], TIs[end]],
	    width=600,
	    height=400,
	    hovermode="x unified",
	)
	fig3
end


# ╔═╡ eb0f1e07-24a1-4aa4-b3ef-b63caa97de34
begin
	# Reducing tissues's signal
	signal_myoc5 = reshape(
	    mean(abs.(mean(mag5bb[:, spins[1], :, :], 3)), 1), length(T2ps), length(RRs)
	)
	signal_bloo5 = reshape(
	    mean(abs.(mean(mag5bb[:, spins[2], :, :], 3)), 1), length(T2ps), length(RRs)
	)
	diff_bloo_myoc5 = abs.(signal_bloo5 .- signal_myoc5)
	# Mean
	mean_myoc5 = mean(signal_myoc5, 2)
	mean_bloo5 = mean(signal_bloo5, 2)
	mean_diff5 = mean(diff_bloo_myoc5,2)
	# Std
	std_myoc5  = std(signal_myoc5, 2)
	std_bloo5  = std(signal_bloo5, 2)
	std_diff5 = std(diff_bloo_myoc5,2)
	# Plotting results
	# Mean
	s15 = scatter(;
	    x=T2ps,
	    y=mean_myoc5[:],
	    name=labels[1],
	    legendgroup=labels[1],
	    line=attr(; color=colors[1]),
	)
	s25 = scatter(;
	    x=T2ps,
	    y=mean_bloo5[:],
	    name=labels[2],
	    legendgroup=labels[2],
	    line=attr(; color=colors[2]),
	)
	s35 = scatter(;
		x=T2ps,
		y=mean_diff5[:],
		name="|Blood-Myoc|",
		legendgroup="|Blood-Myoc|",
		line=attr(color=colors[4])
	)
	# Std
	s45 = scatter(;
	    x=[T2ps; reverse(T2ps)],
	    y=[(mean_myoc5 .- std_myoc5)[:]; reverse((mean_myoc5 .+ std_myoc5)[:])],
	    name=labels[1],
	    legendgroup=labels[1],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,0,255,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s55 = scatter(;
	    x=[T2ps; reverse(T2ps)],
	    y=[(mean_bloo5 .- std_bloo5)[:]; reverse((mean_bloo5 .+ std_bloo5)[:])],
	    name=labels[2],
	    legendgroup=labels[2],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(255,0,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s65 = scatter(;
		x=[T2ps; reverse(T2ps)],
		y=[(mean_diff5 .- std_diff5)[:]; reverse((mean_diff5 .+ std_diff5)[:])],
		name="|Blood-Myoc|",legendgroup="|Blood-Myoc|",
		showlegend=false,
		fill="toself",
		fillcolor="rgba(255,0,255,0.2)",
		line=attr(color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	# Plots
	fig5 = plot([s15, s25, s35, s45, s55, s65])
	relayout!(
	    fig5;
	    yaxis=attr(; title="Signal [a.u.]", tickmode="array"),
	    xaxis=attr(;
	        title="T2prep duration [ms]",
	        tickmode="array",
	        tickvals=[T2ps[1], 50, 70, T2ps[end]],
	        constrain="domain",
	    ),
	    font=attr(; family="CMU Serif", size=16, scaleanchor="x", scaleratio=1),
	    yaxis_range=[0, 0.2],
		xaxis_range=[T2ps[1], T2ps[end]],
	    width=600,
	    height=400,
	    hovermode="x unified",
	)
	fig5
end

# ╔═╡ 10e68e25-9f68-46c0-b1b0-e9b124706b67
begin
	# Reducing tissues's signal
	signal_myocbb = reshape(
	    mean(abs.(mean(mag5bb[:, spins[1], :, :], 3)), 1), length(T2ps), length(RRs)
	)
	signal_bloodbb = reshape(
	    mean(abs.(mean(mag5bb[:, spins[2], :, :], 3)), 1), length(T2ps), length(RRs)
	)
	signal_myocrf = reshape(
	    mean(abs.(mean(mag5rf[:, spins[1], :, :], 3)), 1), length(T2ps), length(RRs)
	)
	signal_bloodrf = reshape(
	    mean(abs.(mean(mag5rf[:, spins[2], :, :], 3)), 1), length(T2ps), length(RRs)
	)
	# Substracted
	signal_myoc6 = abs.(signal_myocrf .- signal_myocbb)
	signal_bloo6 = abs.(signal_bloodrf .- signal_bloodbb)
	diff_bloo_myoc6 = abs.(signal_bloo6 .- signal_myoc6)
	# Mean
	mean_myoc6 = mean(signal_myoc6, 2)
	mean_bloo6 = mean(signal_bloo6, 2)
	mean_diff6 = mean(diff_bloo_myoc6,2)
	# Std
	std_myoc6  = std(signal_myoc6, 2)
	std_bloo6  = std(signal_bloo6, 2)
	std_diff6 = std(diff_bloo_myoc6,2)
	# Plotting results
	# Mean
	s16 = scatter(;
	    x=T2ps,
	    y=mean_myoc6[:],
	    name="Substracted "*labels[1],
	    legendgroup=labels[1],
	    line=attr(; color=colors[1]),
	)
	s26 = scatter(;
	    x=T2ps,
	    y=mean_bloo6[:],
	    name="Substracted "*labels[2],
	    legendgroup=labels[2],
	    line=attr(; color=colors[2]),
	)
	s36 = scatter(;
		x=T2ps,
		y=mean_diff6[:],
		name="|Sub. Blood - Sub. Myoc|",
		legendgroup="|Blood-Myoc|",
		line=attr(color=colors[4])
	)
	# Std
	s46 = scatter(;
	    x=[T2ps; reverse(T2ps)],
	    y=[(mean_myoc6 .- std_myoc6)[:]; reverse((mean_myoc6 .+ std_myoc6)[:])],
	    name="Substracted "*labels[1],
	    legendgroup=labels[1],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(0,0,255,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s56 = scatter(;
	    x=[T2ps; reverse(T2ps)],
	    y=[(mean_bloo6 .- std_bloo6)[:]; reverse((mean_bloo6 .+ std_bloo6)[:])],
	    name="Substracted "*labels[2],
	    legendgroup=labels[2],
	    showlegend=false,
	    fill="toself",
	    fillcolor="rgba(255,0,0,0.2)",
	    line=attr(; color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	s66 = scatter(;
		x=[T2ps; reverse(T2ps)],
		y=[(mean_diff6 .- std_diff6)[:]; reverse((mean_diff6 .+ std_diff6)[:])],
		name="|Sub. Blood - Sub. Myoc|",legendgroup="|Blood-Myoc|",
		showlegend=false,
		fill="toself",
		fillcolor="rgba(255,0,255,0.2)",
		line=attr(color="rgba(0,0,0,0)"),
		hoverinfo="none"
	)
	# Plots
	fig6 = plot([s16, s26, s36, s46, s56, s66])
	relayout!(
	    fig6;
	    yaxis=attr(; title="Signal [a.u.]", tickmode="array"),
	    xaxis=attr(;
	        title="T2prep duration [deg]",
	        tickmode="array",
	        tickvals=[T2ps[1], 50, 70, T2ps[end]],
	        constrain="domain",
	    ),
	    font=attr(; family="CMU Serif", size=16, scaleanchor="x", scaleratio=1),
	    yaxis_range=[0, 0.2],
		xaxis_range=[T2ps[1], T2ps[end]],
	    width=600,
	    height=400,
	    hovermode="x unified",
	)
	fig6
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
KomaMRICore = "4baa4f4d-2ae9-40db-8331-a7d1080e3f4e"
KomaMRIPlots = "76db0263-63f3-4d26-bb9a-5dba378db904"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Suppressor = "fd094767-a336-5f1f-9728-57cf17d0bbfb"

[compat]
PlutoUI = "~0.7.58"
ProgressLogging = "~0.1.4"
Suppressor = "~0.2.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "9b71b5fa18cd815278316446dec81dab9113a167"

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
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
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

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires"]
git-tree-sha1 = "e2144b631226d9eeab2d746ca8880b7ccff504ae"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.3"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "cb25e4b105cc927052c2314f8291854ea59bf70a"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.4"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Blink]]
deps = ["Base64", "Distributed", "HTTP", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Pkg", "Reexport", "Sockets", "WebIO"]
git-tree-sha1 = "bc93511973d1f949d45b0ea17878e6cb0ad484a1"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.9"

[[deps.BufferedStreams]]
git-tree-sha1 = "6863c5b7fc997eadcabdbaf6c5f201dc30032643"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.2.2"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

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

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

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

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "64d8e93700c7a3f28f717d265382d52fac9fa1c1"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.12"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dd3b49277ec2bb2c6b94eb1604d4d0616016f7a6"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.11.2+0"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

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

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2ef87eeaa28713cb010f9fb0be288b6c1a4ecd53"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.1.0+0"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "MacroTools", "PrecompileTools", "Requires", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "5126765c5847f74758c411c994312052eb7117ef"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.27"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.KomaMRIBase]]
deps = ["Interpolations", "MAT", "MRIBase", "Parameters", "Reexport"]
git-tree-sha1 = "640eaa85cedeeabbc4809295f27693ab414f9530"
uuid = "d0bc0b20-b151-4d03-b2a4-6ca51751cb9c"
version = "0.9.0"

[[deps.KomaMRICore]]
deps = ["Adapt", "Functors", "KernelAbstractions", "KomaMRIBase", "ProgressMeter", "Reexport", "ThreadsX"]
git-tree-sha1 = "728828e5110ce74a96a0f9780a0ba46640357dfc"
uuid = "4baa4f4d-2ae9-40db-8331-a7d1080e3f4e"
version = "0.9.0"

    [deps.KomaMRICore.extensions]
    KomaAMDGPUExt = "AMDGPU"
    KomaCUDAExt = "CUDA"
    KomaMetalExt = "Metal"
    KomaoneAPIExt = "oneAPI"

    [deps.KomaMRICore.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    oneAPI = "8f75cd03-7ff8-4ecb-9b8f-daf728133b1b"

[[deps.KomaMRIPlots]]
deps = ["Interpolations", "Kaleido_jll", "KomaMRIBase", "MAT", "PlotlyJS", "Reexport"]
git-tree-sha1 = "83270cba62a2660dc49e3f7e704418b3d57f3125"
uuid = "76db0263-63f3-4d26-bb9a-5dba378db904"
version = "0.9.0"
weakdeps = ["PlutoPlotly"]

    [deps.KomaMRIPlots.extensions]
    KomaPlotsPlutoPlotlyExt = "PlutoPlotly"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Requires", "Unicode"]
git-tree-sha1 = "4ad43cb0a4bb5e5b1506e1d1f48646d7e0c80363"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "9.1.2"

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

    [deps.LLVM.weakdeps]
    BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "05a8bd5a42309a9ec82f700876903abce1017dd3"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.34+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "1d2dd9b186742b0f317f2530ddcbf00eebb18e96"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.7"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "19d4bd098928a3263693991500d05d74dbdc2004"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.2.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "fde81c9f9c94fe5fbeaed7b3f1330305cf9a327c"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.0+0"

[[deps.MRIBase]]
deps = ["AbstractNFFTs", "LinearAlgebra", "NFFTTools"]
git-tree-sha1 = "57979500dbdd130fc92359f1ecd6714051ed78eb"
uuid = "f7771a9a-6e57-4e71-863b-6e4b6a2f17df"
version = "0.4.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

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

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f12a29c4400ba812841c6ace3f4efbb6dbb3ba01"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "3b2db451a872b20519ebb0cec759d3d81a1c6bcb"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.20"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "MbedTLS", "Pkg", "Sockets"]
git-tree-sha1 = "7295d849103ac4fcbe3b2e439f229c5cc77b9b69"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "1.0.2"

[[deps.NFFTTools]]
deps = ["AbstractFFTs", "AbstractNFFTs", "FFTW", "LinearAlgebra"]
git-tree-sha1 = "d6a68b7ffbd50b4c99e514a1a6fb8ce84f6e247e"
uuid = "7424e34d-94f7-41d6-98a0-85abaf1b6c91"
version = "0.2.6"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e25c1778a98e34219a00455d6e4384e017ea9762"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.6+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

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
git-tree-sha1 = "f198c8a80c08987a2915156e6e6131e5d73b97f4"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.14"

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
git-tree-sha1 = "3210de4d88af7ca5de9e26305758a59aabc48aac"
uuid = "f2990250-8cf9-495f-b13a-cce12b45703c"
version = "2.2.5"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "653b48f9c4170343c43c2ea0267e451b68d69051"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.5.0"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

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
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "02d31ad62838181c1a3a5fd23a1ce5914a643601"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

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
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Suppressor]]
deps = ["Logging"]
git-tree-sha1 = "6dbb5b635c5437c68c28c2ac9e39b87138f37c0a"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.8"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

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

[[deps.TranscodingStreams]]
git-tree-sha1 = "e84b3a11b9bece70d14cce63406bbc79ed3464d2"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.2"

[[deps.Transducers]]
deps = ["Accessors", "Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "SplittablesBase", "Tables"]
git-tree-sha1 = "5215a069867476fc8e3469602006b9670e68da23"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.82"

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
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

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

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "2d17fabcd17e67d7625ce9c531fb9f40b7c42ce4"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.2.1"

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
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

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
# ╟─42417eba-46e0-400e-a874-299533362c41
# ╠═0b6c1f72-b040-483c-969b-88bfe09b32c3
# ╟─4aab8a6f-d2ba-46f0-a1bf-ddc0dcf8437f
# ╠═88eb41a5-d8c2-4f0e-b379-a8b05a341a82
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
# ╟─b0b53632-e6b8-47e9-8d3e-79f1e599315a
# ╠═14cf6859-a4b9-4671-b022-659781c55144
# ╟─6649f46e-3565-4cd6-86a5-9b03d00cc3db
# ╟─c64d38cb-2433-4b90-abbf-1fc938f70584
# ╠═c952ecf1-25ef-4b48-9c8f-e53ded302629
# ╟─a8751a97-c131-407b-b211-573660452142
# ╟─85834365-238b-4193-a0c0-d1859416ab6c
# ╠═50b551d8-a212-4319-b25b-1f8ecdd2a491
# ╟─eb0f1e07-24a1-4aa4-b3ef-b63caa97de34
# ╟─4b1dc8bf-feb1-42b4-94f7-1b27975d944c
# ╟─10e68e25-9f68-46c0-b1b0-e9b124706b67
# ╟─8dc11175-ebc8-407a-ab0b-6d543f849a72
# ╟─3d7e7d20-a77a-48b3-ad2e-6b621227be16
# ╟─abea2c43-d83e-4438-8cd3-4be06b8174b3
# ╟─7deadd58-b202-4508-b4c7-686f742cb713
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

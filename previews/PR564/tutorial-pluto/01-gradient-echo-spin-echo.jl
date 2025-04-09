### A Pluto.jl notebook ###
# v0.20.5

#> [frontmatter]
#> title = "Understanding basic MRI sequences"
#> tags = ["educational"]
#> description = "Free Induction Decay (FID), Gradient Echo (GE), and Spin Echo (SE)"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

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

# ╔═╡ d6b1729a-874d-11ee-151a-9b0fcce2c4fd
using KomaMRICore, KomaMRIPlots, FFTW, PlutoPlotly, PlutoUI

# ╔═╡ 5df97874-f09c-4173-a2f6-893db322ccaf
md"# Understanding basic MRI sequences"

# ╔═╡ 8529f36d-2d39-4b45-a821-01c8346539fd
TableOfContents() # There should be a table of contents on the right --->

# ╔═╡ 6dfe338d-de85-4adb-b030-09455fae78a0
md"""
Welcome to the hands-on session on MRI simulation. Let's have some fun!

If you have any doubts about how to use a function, please search in the **Live Docs** at the bottom right.
"""

# ╔═╡ 8e474add-8651-431b-b481-7a139037dbd2
md"""# 1. Free Induction Decay (FID)
The free induction decay is the simplest observable NMR signal. This signal is the one that follows a single tipping RF pulse.

$(PlutoUI.Resource("https://raw.githubusercontent.com/LIBREhub/MRI-processing-2023/main/02-simulation/Figures/FID.png", :width=>"300px"))
To recreate this experiment, we will need to define a `Sequence`:
 - (1.1) A 90-deg block RF pulse, put it in a variable `seq` (check `PulseDesigner.RF_hard`'s docs using the Live Docs)
 - (1.2) An ADC to capture the signal in a variable `adc`, concatenate with (1.1) using `seq += adc`
 - (1.3) Plot the generated `Sequence` (check `plot_seq`'s docs)

For the hardware limits use the default scanner `sys = Scanner()`. 
"""

# ╔═╡ c6e33cb8-f42c-4643-9257-124d2804d3da
# (1.1) A 90-deg block RF pulse
begin
	sys = Scanner()
	durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse
	rf = PulseDesigner.RF_hard(sys.B1, durRF, sys)
end

# ╔═╡ 0266632d-5ca4-4196-a523-33a66dd70e0c
# (1.2) An ADC to capture the signal
adc = ADC(100, 50e-3)

# ╔═╡ 0975547d-67d9-4e6b-88ff-a9dd06a7f9ef
# (1.3) Plot the generated Sequence
begin
	seq = Sequence()
	seq += rf
	seq += adc
	plot_seq(seq; slider=false)
end

# ╔═╡ f11a2fa2-eff9-4979-b739-3da2b24a9a45
md"""
Generate a virtual object:
 - (1.4) A Phantom with 20 spins, with properties:
   - `obj.x`  = [-1, 1] mm (20 spins along the $x$-axis)
   - `obj.ρ`  .= 1
   - `obj.T1` .= 500 ms
   - `obj.T2` .= 50 ms
 - (1.5) Plot the generated `Phantom` (check `plot_phantom_map`'s docs)
"""

# ╔═╡ d16efa62-dce7-4ec3-9e3c-b5e1677377fc
# (1.4) A Phantom with 20 spins
begin
	obj = Phantom(x=collect(range(-1e-3,1e-3,20)))
	obj.ρ .= 1
	obj.T1 .= 500e-3
	obj.T2 .= 50e-3
	nothing
end

# ╔═╡ 35ff3402-dc36-4b91-bec9-b4d21faf3e68
# (1.5) Plot the generated Phantom
plot_phantom_map(obj, :T1)

# ╔═╡ ea542271-01c2-4962-a708-804b23a861b9
md"""
 - (1.6) Finally, use the generated `seq`, `obj`, and `sys` to simulate the FID (check `simulate`'s docs)
 - (1.7) Plot the resulting raw data with `plot_signal`.
 - (1.8) Is the signal the same as `plot(t, exp.(-t ./ T2))`?
"""

# ╔═╡ c47a50b8-c930-4c96-9b34-2772186634d9
# (1.6) Finally, use the generated seq, obj, and sys to simulate the FID
raw = simulate(obj, seq, sys)

# ╔═╡ 7a66ab47-918f-4582-895f-1b4690562051
# (1.7) Plot the resulting raw data with plot_signal
plot_signal(raw; slider=false)

# ╔═╡ 1231b832-47b1-4ccb-9b56-a67838598cc7
# (1.8) Is the signal the same as `plot(t, exp.(-t ./ T2))`?
begin
	t = range(0, 50, 100)
	t2_decay(t) = scatter(x=t, y=20.0.*exp.(-t ./ 50), name="T2-decay", marker_color="purple")
	plot(t2_decay(t), Layout(yaxis_range=[0, 20.1]))
end

# ╔═╡ e4c80c24-20fd-42e5-9dcd-a65958569c01
md"""
# 2. Gradient Echo

$(Resource("https://raw.githubusercontent.com/LIBREhub/MRI-processing-2023/main/02-simulation/Figures/GRE.gif", :width=>"400px"))

The gradient echo is one of the first steps to create an image. The big 
breakthrough was the addition of linearly increasing magnetic fields, or gradients, to encode the spin's positions in their frequency (Mmmh, someone said Fourier?). This works due to the fact that the frequency $$f$$ of a spin is 

$$f(x) = \frac{\gamma}{2\pi} B_z(x) = \frac{\gamma}{2\pi} G_x x.$$

Let's create a different sequence.
 - Create a 90-deg hard RF pulse and put it in a variable `seq_gre`
 - (2.1) Create a gradient with area `-Ax` using `gx_pre = Grad(A,T,rise,fall)` append to `seq_gre`. As an optional challenge, put `gx_pre.rise` and `gx_pre.fall` so the satisfy the `sys` requierements
 - (2.2) Append a `Sequence` block called `readout` that includes: 
   - A gradient of twice the area, or `2Ax`. Call it `gx`
   - An `ADC` with `adc2.delay = gx.rise` and `adc2.T = gx.T`
 - (2.3) Plot `seq_gre` and its k-space
 - (2.4) Plot the $$k$$-space with the `plot_kspace` function
"""

# ╔═╡ 74666c1a-2673-4936-982b-6229bf92af66
md"""
 - (2.5) Simulate the `seq_gre` sequence
 - (2.6) Plot the simulated signal
 - (2.7) Reconstruct the 1D image
 - (2.8) Do you notice anything weird? If the answer is yes, try adjusting `Ax` to change the `FOV` of the acquisition
"""

# ╔═╡ 0f96a83d-96ef-4768-9330-87c466e35c93
# (2.8) Do you notice anything weird? Change Ax!
@bind Ax Slider(range(0, 20, 20)*1e-5, default=10e-5) # Gradient's area in [T/m s]

# ╔═╡ 9179aa40-bb40-4a36-ae1e-00ae42935a5f
# (2.1) Create a gradient `gx_pre`, use the variable `Ax`!!
begin
	T_gx_pre = 10e-3
	gx_pre = Grad(-Ax/T_gx_pre, T_gx_pre, 0, 0)
	seq_gre = Sequence()
	seq_gre += rf
	seq_gre += gx_pre
# (2.2) Append a `Sequence` block called `readout`
	gx = Grad(2*Ax/(2T_gx_pre), 2T_gx_pre, 0, 0)
	adc2 = ADC(100, 2T_gx_pre)
	readout = Sequence([gx;;], [RF(0,0);;], [adc2])
	seq_gre += readout
end

# ╔═╡ 8b4a1ad9-2d6a-4c8f-bb8e-f43c2d058195
# (2.3) Plot `seq_gre` and the k-space
plot_seq(seq_gre; slider=false)

# ╔═╡ 3abca406-2e6b-4b37-8835-65cfad9d0caa
# (2.4) Plot the $k$-space with the `plot_kspace` function
kspace_gre = plot_kspace(seq_gre)

# ╔═╡ ada602d2-4f4b-4fb4-a763-8a639e05ff38
# (2.5) Simulate the seq_gre sequence
raw_gre = simulate(obj, seq_gre, sys)

# ╔═╡ 9a88a54b-bcc7-41ad-8e60-f4d450dccb2d
# (2.7) Reconstruct the 1D image
begin
    fftc(x; dims=[1,2]) = fftshift(fft(ifftshift(x, dims), dims), dims)/prod(size(x)[dims])
    recon_gre = plot(abs.(fftc(raw_gre.profiles[1].data)))
end

# ╔═╡ 41d14dec-b852-4316-aefb-c3d08fa43216
# (2.6) Plot the simulated signal
begin
	t_adc_gre = KomaMRICore.get_adc_sampling_times(seq_gre)*1e3
	signal_gre = plot_signal(raw_gre; slider=false)
    addtraces!(signal_gre, t2_decay(t_adc_gre))
	signal_gre
end

# ╔═╡ 97104c46-e81f-444a-957f-0bbb1b02f1b8
md"""
# 3. $T_{2}^{*}$-decay

The $$T_{2}^{*}$$-decay is the signal decay produced by microscopic distribution of off-resonance.

$(Resource("https://raw.githubusercontent.com/LIBREhub/MRI-processing-2023/main/02-simulation/Figures/T2star.png", :width=>"400px"))

The exact distribution of off-resonance is

$$p_{\Delta w}(w) = \frac{T_2^{'}}{\pi(1+T_2^{'2} w^2)},\quad\text{with }\frac{1}{T_2^{*}} = \frac{1}{T_2} + \frac{1}{T_2^{'}}.$$

In this excercise we will simplify this distribution, but we will obtain a similar effect.

- (3.1) Create a new phantom named `obj_t2star` with spins at the same positions as the original phantom `obj`, each having a linear distribution of off-resonance. To achieve this, follow these steps:
   * (3.1.1) Create an empty phantom called `obj_t2star`.
   * (3.1.2) Create a linear off-resonance distribution such that the range $$2\pi [-10, 10]\,\mathrm{rad/s}$$ is covered uniformly with $$N_{\mathrm{isochromats}} = 20$$ (use the function `range(start, stop, length)`).
   * (3.1.3) Iterate over the elements `off` of the linear distribution (`for` loop) and create copies of the original phantom (`obj_aux = copy(obj)`) and set the off-resonance of that copy to `off` with `obj_aux.Δw .= off`.
   * (3.1.4) Update `obj_t2star` by appending the modified copies `obj_aux` (`obj_t2star += obj_aux`).
   * (3.1.5) Finally, outside the loop, divide the proton density `obj_t2star.ρ` by $$N_{\mathrm{isochromats}} = 20$$ and rename the phantom `obj_t2star.name = "T2 star phantom"`.

 - (3.2) Plot `obj_t2star` with `plot_phantom_map(obj_t2star, :Δw)` and verify it is correct

"""

# ╔═╡ ee7e81e7-484c-44a8-a191-f73e24707ce9
# (3.1) Create the new obj_t2star phantom 
begin
    # (3.1.1) Create an empty phantom
	obj_t2star = Phantom()
    # (3.1.2) Define the linear off-resonance distribution
	Niso = 20
	linear_offresonance_distribution = 2π .* range(-10, 10, Niso)
	# (3.1.3) Iterate over the linear off-resonance distribution and ...
	for off = linear_offresonance_distribution
		# ... copy the original phantom and modify its off-resonance
	    aux = copy(obj)
		aux.Δw .= off
		aux.y  .+= off * 1e-6  # So the distribution is visible
		# (3.1.4) Update the phantom
		obj_t2star += aux
	end
	# (3.1.5) Divide the proton density and rename the phantom
	obj_t2star.ρ .= 1.0 / Niso
	obj_t2star.name = "T2 star phantom"
end

# ╔═╡ 2ee7ba47-02e5-4b02-a162-ddbd5ed47c7b
# (3.2) Plot obj_t2star
plot_phantom_map(obj_t2star, :Δw)

# ╔═╡ 27686262-1a1e-45fa-b4ee-90ae1d9ee34e
md"""
 - (3.3) Simulate the `seq_gre` sequence
 - (3.4) Plot the simulated signal
 - (3.5) Compare the plot in (3.5) with (2.6)
 - (3.6) Reconstruct the 1D image
"""

# ╔═╡ e4ef5145-a63c-4f91-ac04-3b5bf16c0842
# (3.3) Simulate the seq_gre sequence
raw_t2_star_gre = simulate(obj_t2star, seq_gre, sys)

# ╔═╡ 1a83d897-705b-443d-89a4-ea5e3e6a3c07
# (3.4) Plot the simulated signal
begin
	signal_t2_star_gre = plot_signal(raw_t2_star_gre; slider=false)
	addtraces!(signal_t2_star_gre, t2_decay(t_adc_gre))
	signal_t2_star_gre
end

# ╔═╡ 18c82ff1-0bde-4fa0-848c-d0eb73d1ac7c
# (3.5) Compare the plot in (3.4) with (2.6)
begin
	signal_layout = Layout(yaxis=attr(range=[-5, 16.5]))
	relayout!(signal_gre, signal_layout; title="GRE-T2")
	relayout!(signal_t2_star_gre, signal_layout; title="GRE-T2*")
	fig_signal_2 = [signal_gre signal_t2_star_gre]
	relayout(fig_signal_2, showlegend=false)
end

# ╔═╡ 4a4a6bd3-b820-479c-89e3-f3ce79a316db
# (3.6) Reconstruct the 1D image
recon_t2_star_gre = plot(abs.(fftc(raw_t2_star_gre.profiles[1].data)))

# ╔═╡ 964404f6-7f46-4df9-ad98-921948c3be69
begin
	recon_layout = Layout(yaxis=attr(range=[0, 0.8]))
	relayout!(recon_gre, recon_layout; title="GRE-T2")
	relayout!(recon_t2_star_gre, recon_layout; title="GRE-T2*")
	fig_recon_2 = [recon_gre recon_t2_star_gre]
	relayout(fig_recon_2, showlegend=false)
end

# ╔═╡ 3357a283-a234-4d15-8fdf-7fbec58b33a7
md"""
# 4. Spin Echo

$(Resource("https://raw.githubusercontent.com/LIBREhub/MRI-processing-2023/main/02-simulation/Figures/SE.gif", :width=>"400px"))

The spin echo experiment has the advantage that the echo signal amplitud it is modulated by $$\exp(-t/T_2)$$ and not $$\exp(-t/T_2^{*})$$.

For this section we will use the phantom `obj_t2star` and a new sequence `seq_se`.

For this sequence we will need:
 - (4.1) A 90deg hard RF pulse
 - (4.2) A `Delay` of $$\mathrm{TE}/2$$ with a positive gradient (area `Ax`)
 - (4.3) A 180deg hard RF pulse
 - (4.4) A readout gradient of area `2Ax` with an ADC (similar to (2.2)), such that the middle of the gradient and ADC are in $$\mathrm{TE}$$
 - (4.5) Create concatenating these blocks into a sequence called `seq_se`
 - (4.6) Plot `seq_se` and its k-space. Is the k-space the same as `seq_gre` in (2.3)?
"""

# ╔═╡ 27e65680-22a0-4079-b6df-d60a3218e52e
# (4.5) Create concatenating these blocks into a sequence called `seq_se`
begin
	# (4.1) A 90deg hard RF pulse
	seq_se = Sequence()
	seq_se += rf
    # (4.2) A `Delay` of TE/2 with a positive gradient (area `Ax`)
	seq_se += -1*gx_pre
	# (4.3) A 180deg hard RF pulse
	seq_se += (0.0+2.0im)*rf
	# (4.4) A readout gradient of area `2Ax` with an ADC (similar to (2.2)), such that the middle of the gradient and ADC are in $$\mathrm{TE}$$
	seq_se += readout
end

# ╔═╡ f1f3b700-5916-496f-b938-46f7f08b4eb6
# (4.6) Plot seq_se and its k-space. Is the k-space the same as seq_gre in (2.3)?
plot_seq(seq_se; slider=false)

# ╔═╡ 4e1434e1-673f-4206-a271-9edec10ebd6a
kspace_se = plot_kspace(seq_se)

# ╔═╡ c02f3898-10cb-4f1e-b5ef-eb42b803baed
begin
	relayout!(kspace_gre; title="GRE")
	relayout!(kspace_se; title="SE")
	fig_kspace = [kspace_gre kspace_se]
	relayout(fig_kspace, showlegend=false)
end

# ╔═╡ 45952512-aaf1-43d8-a95e-c32bb2633f42
md"""
 - (4.7) Simulate using `seq_se` and `obj_t2star`
 - (4.8) Compare the signal obtained in (4.6) with the one at (3.5)
 - (4.9) Reconstruct the 1D image
"""

# ╔═╡ 97479437-9ce3-4b33-9134-0f2af89bccb5
# (4.7) Simulate using seq_se and obj_t2star
raw_t2_star_se = simulate(obj_t2star, seq_se, sys)

# ╔═╡ 1c79b37e-d4e0-490f-9466-20ce28f017ae
# (4.8) Compare the signal obtained in (4.7) with the one at (3.4)
begin
	t_adc_se = KomaMRICore.get_adc_sampling_times(seq_se)*1e3
	signal_t2_star_se = plot_signal(raw_t2_star_se; slider=false)
	addtraces!(signal_t2_star_se, t2_decay(t_adc_se))
	relayout!(signal_t2_star_se, signal_layout; title="SE")
	fig_signal_3 = [signal_gre signal_t2_star_gre signal_t2_star_se]
	relayout(fig_signal_3, showlegend=false)
end

# ╔═╡ 2e65ae31-f50a-462b-9744-80bf6cdb388e
# (4.9) Reconstruct the 1D image
recon_t2_star_se = plot(abs.(fftc(raw_t2_star_se.profiles[1].data)))

# ╔═╡ 34824db7-13c4-45e2-befa-f027b9b585c0
begin
	relayout!(recon_t2_star_se, recon_layout; title="SE")
	fig_recon_3 = [recon_gre recon_t2_star_gre recon_t2_star_se]
	relayout(fig_recon_3, showlegend=false)
end

# ╔═╡ fe8bbcd2-e8f5-4225-80c3-47e73176fb3d
md"""
Congratulations! you finished the simulation hands-on session 🥳!
"""

# ╔═╡ ab8dc1ce-d1ef-43a0-9495-dac931b52aec
# Set this boolean to `true` when you finish
activity_finished = true

# ╔═╡ 58be4150-2b7a-4f9e-a7d7-40a086fd3a53
if activity_finished
    html"""
    <script>
    const {default: confetti} = await import("https://cdn.skypack.dev/canvas-confetti@1")
    confetti()
    </script>
    """
end

# ╔═╡ abea2c43-d83e-4438-8cd3-4be06b8174b3
md"""# Reproducibility

This [Pluto notebook](https://plutojl.org/) is reproducible by default, as it has an embedded `Project.toml` and `Manifest.toml`, that store the exact package versions used to create the notebook."""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
KomaMRICore = "4baa4f4d-2ae9-40db-8331-a7d1080e3f4e"
KomaMRIPlots = "76db0263-63f3-4d26-bb9a-5dba378db904"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FFTW = "~1.8.1"
KomaMRICore = "~0.9.1"
KomaMRIPlots = "~0.9.2"
PlutoPlotly = "~0.5.0"
PlutoUI = "~0.7.62"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "02f4cbcff5318906be68878937be9de74cc2ea22"

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
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgCheck]]
git-tree-sha1 = "f9e9a66c9b7be1ad7372bbd9b062d9230c30c5ce"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "b5bb4dc6248fde467be2a863eb8452993e74d402"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "1.1.1"

    [deps.Atomix.extensions]
    AtomixCUDAExt = "CUDA"
    AtomixMetalExt = "Metal"
    AtomixOpenCLExt = "OpenCL"
    AtomixoneAPIExt = "oneAPI"

    [deps.Atomix.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    OpenCL = "08131aa3-fb12-5dee-8b74-c09406e224a2"
    oneAPI = "8f75cd03-7ff8-4ecb-9b8f-daf728133b1b"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra"]
git-tree-sha1 = "26f41e1df02c330c4fa1e98d4aa2168fdafc9b1f"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.4"

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
version = "1.11.0"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

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
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

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
version = "1.11.0"

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
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

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
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

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
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "c67b33b085f6e2faf8bf79a61962e7339a81129c"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.15"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f93a9ce66cd89c9ba7a4695a47fd93b4c6bc59fa"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.0+0"

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
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

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
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

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
deps = ["Adapt", "Atomix", "InteractiveUtils", "MacroTools", "PrecompileTools", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "80d268b2f4e396edc5ea004d1e0f569231c71e9e"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.34"

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
git-tree-sha1 = "636d875ddea2355c42901d67e345069618aa8b01"
uuid = "d0bc0b20-b151-4d03-b2a4-6ca51751cb9c"
version = "0.9.1"

[[deps.KomaMRICore]]
deps = ["Adapt", "Functors", "KernelAbstractions", "KomaMRIBase", "ProgressMeter", "Reexport", "ThreadsX"]
git-tree-sha1 = "08af8891306737ae3abbe0cc333d553220593392"
uuid = "4baa4f4d-2ae9-40db-8331-a7d1080e3f4e"
version = "0.9.1"

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
deps = ["Interpolations", "Kaleido_jll", "KomaMRIBase", "MAT", "PlotlyJS", "QMRIColors", "Reexport"]
git-tree-sha1 = "a5d24d52b960d2a6eefd15abe038b881ea4c21ff"
uuid = "76db0263-63f3-4d26-bb9a-5dba378db904"
version = "0.9.2"
weakdeps = ["PlutoPlotly"]

    [deps.KomaMRIPlots.extensions]
    KomaPlotsPlutoPlotlyExt = "PlutoPlotly"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MAT]]
deps = ["BufferedStreams", "CodecZlib", "HDF5", "SparseArrays"]
git-tree-sha1 = "1d2dd9b186742b0f317f2530ddcbf00eebb18e96"
uuid = "23992714-dd62-5051-b70f-ba57cb901cac"
version = "0.10.7"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "3aa3210044138a1749dbd350a9ba8680869eb503"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.0+1"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "ff91ca13c7c472cef700f301c8d752bc2aaff1a8"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.3+0"

[[deps.MRIBase]]
deps = ["AbstractNFFTs", "LinearAlgebra", "NFFTTools"]
git-tree-sha1 = "57979500dbdd130fc92359f1ecd6714051ed78eb"
uuid = "f7771a9a-6e57-4e71-863b-6e4b6a2f17df"
version = "0.4.4"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

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
git-tree-sha1 = "a414039192a155fb38c4599a60110f0018c6ec82"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.16.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "047b66eb62f3cae59ed260ebb9075a32a04350f1"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.7+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a9697f1d06cc3eb3fb3ad49cc67f2cfabaac31ea"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.16+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

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
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "90af5c9238c1b3b25421f1fdfffd1e8fca7a7133"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.20"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "PlotlyKaleido", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "e415b25fdec06e57590a7d5ac8e0cf662fa317e2"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.15"

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
deps = ["Artifacts", "Base64", "JSON", "Kaleido_jll"]
git-tree-sha1 = "9ef5c9e588ec7e912f01a76c7fd3dddf1913d4f2"
uuid = "f2990250-8cf9-495f-b13a-cce12b45703c"
version = "2.3.0"

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
git-tree-sha1 = "d3de2694b52a01ce61a036f18ea9c0f61c4a9230"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.62"

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
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.QMRIColors]]
deps = ["Colors", "DelimitedFiles"]
git-tree-sha1 = "77d0d5f446d544beaa568b837fade390d9d38721"
uuid = "2bec176e-9e8c-4764-a62f-295118d1ec05"
version = "1.0.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

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
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

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
version = "1.11.0"

[[deps.ThreadsX]]
deps = ["Accessors", "ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "SplittablesBase", "Transducers"]
git-tree-sha1 = "70bd8244f4834d46c3d68bd09e7792d8f571ef04"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.12"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Transducers]]
deps = ["Accessors", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "SplittablesBase", "Tables"]
git-tree-sha1 = "7deeab4ff96b85c5f72c824cae53a1398da3d1cb"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.84"

    [deps.Transducers.extensions]
    TransducersAdaptExt = "Adapt"
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "b13c4edda90890e5b04ba24e20a310fbe6f249ff"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.3.0"

    [deps.UnsafeAtomics.extensions]
    UnsafeAtomicsLLVM = ["LLVM"]

    [deps.UnsafeAtomics.weakdeps]
    LLVM = "929cbde3-209d-540e-8aea-75f648917ca0"

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
git-tree-sha1 = "e9aeb174f95385de31e70bd15fa066a505ea82b9"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.7"

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
git-tree-sha1 = "f5733a5a9047722470b95a81e1b172383971105c"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.3+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─5df97874-f09c-4173-a2f6-893db322ccaf
# ╟─8529f36d-2d39-4b45-a821-01c8346539fd
# ╟─6dfe338d-de85-4adb-b030-09455fae78a0
# ╠═d6b1729a-874d-11ee-151a-9b0fcce2c4fd
# ╟─8e474add-8651-431b-b481-7a139037dbd2
# ╠═c6e33cb8-f42c-4643-9257-124d2804d3da
# ╠═0266632d-5ca4-4196-a523-33a66dd70e0c
# ╠═0975547d-67d9-4e6b-88ff-a9dd06a7f9ef
# ╟─f11a2fa2-eff9-4979-b739-3da2b24a9a45
# ╠═d16efa62-dce7-4ec3-9e3c-b5e1677377fc
# ╠═35ff3402-dc36-4b91-bec9-b4d21faf3e68
# ╟─ea542271-01c2-4962-a708-804b23a861b9
# ╠═c47a50b8-c930-4c96-9b34-2772186634d9
# ╠═7a66ab47-918f-4582-895f-1b4690562051
# ╠═1231b832-47b1-4ccb-9b56-a67838598cc7
# ╟─e4c80c24-20fd-42e5-9dcd-a65958569c01
# ╠═9179aa40-bb40-4a36-ae1e-00ae42935a5f
# ╠═8b4a1ad9-2d6a-4c8f-bb8e-f43c2d058195
# ╠═3abca406-2e6b-4b37-8835-65cfad9d0caa
# ╟─74666c1a-2673-4936-982b-6229bf92af66
# ╠═ada602d2-4f4b-4fb4-a763-8a639e05ff38
# ╠═41d14dec-b852-4316-aefb-c3d08fa43216
# ╠═9a88a54b-bcc7-41ad-8e60-f4d450dccb2d
# ╠═0f96a83d-96ef-4768-9330-87c466e35c93
# ╟─97104c46-e81f-444a-957f-0bbb1b02f1b8
# ╠═ee7e81e7-484c-44a8-a191-f73e24707ce9
# ╠═2ee7ba47-02e5-4b02-a162-ddbd5ed47c7b
# ╟─27686262-1a1e-45fa-b4ee-90ae1d9ee34e
# ╠═e4ef5145-a63c-4f91-ac04-3b5bf16c0842
# ╠═1a83d897-705b-443d-89a4-ea5e3e6a3c07
# ╠═18c82ff1-0bde-4fa0-848c-d0eb73d1ac7c
# ╠═4a4a6bd3-b820-479c-89e3-f3ce79a316db
# ╠═964404f6-7f46-4df9-ad98-921948c3be69
# ╟─3357a283-a234-4d15-8fdf-7fbec58b33a7
# ╠═27e65680-22a0-4079-b6df-d60a3218e52e
# ╠═f1f3b700-5916-496f-b938-46f7f08b4eb6
# ╠═4e1434e1-673f-4206-a271-9edec10ebd6a
# ╠═c02f3898-10cb-4f1e-b5ef-eb42b803baed
# ╟─45952512-aaf1-43d8-a95e-c32bb2633f42
# ╠═97479437-9ce3-4b33-9134-0f2af89bccb5
# ╠═1c79b37e-d4e0-490f-9466-20ce28f017ae
# ╠═2e65ae31-f50a-462b-9744-80bf6cdb388e
# ╠═34824db7-13c4-45e2-befa-f027b9b585c0
# ╟─fe8bbcd2-e8f5-4225-80c3-47e73176fb3d
# ╠═ab8dc1ce-d1ef-43a0-9495-dac931b52aec
# ╟─58be4150-2b7a-4f9e-a7d7-40a086fd3a53
# ╟─abea2c43-d83e-4438-8cd3-4be06b8174b3
# ╟─7deadd58-b202-4508-b4c7-686f742cb713
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

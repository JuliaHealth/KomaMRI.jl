#"""
#    phantom = heart_phantom(α=1, β=1, γ=1, fat_bool::Bool=false)
#
#Heart-like LV phantom. The variable `α` is for streching, `β` for contraction, and `γ` for
#rotation.
#
## Arguments
#- `α`: (`::Real`, `=1`) streching parameter
#- `β`: (`::Real`, `=1`) contraction parameter
#- `γ`: (`::Real`, `=1`) rotation parameter
#- `fat_bool`: (`::Bool`, `=false`) fat boolean parameter
#
## Returns
#- `phantom`: (`::Phantom`) Heart-like LV phantom struct
#"""
#heart_phantom(α=1, β=1, γ=1, fat_bool::Bool=false) = begin
#	#PARAMETERS
#	FOV = 10e-2 #m Diameter ventricule
#	N = 21
#	Δxr = FOV/(N-1) #Aprox rec resolution, use Δx_pix and Δy_pix
#	Ns = 50 #number of spins per voxel
#	Δx = Δxr/sqrt(Ns) #spin separation
#	#POSITIONS
#	x = y = -FOV/2:Δx:FOV/2-Δx #spin coordinates
#	x, y = x .+ y'*0, x*0 .+ y' #grid points
#	#PHANTOM
#	⚪(R) =  (x.^2 .+ y.^2 .<= R^2)*1. #Circle of radius R
#	v = FOV/4 #m/s 1/16 th of the FOV during acquisition
#	ωHR = 2π/1 #One heart-beat in one second
#
#	# θ(t) = -π/4*γ*(sin.(ωHR*t).*(sin.(ωHR*t).>0)+0.25.*sin.(ωHR*t).*(sin.(ωHR*t).<0) )
#	ux(x,y,t) = begin
#		strech = 0 #α * v * (x.^2 .+ y.^2) / (FOV/2)^2 .* sign.(x)
#		contract = - β * v * x / (FOV/2)  #expand
#		rotate = - γ * v * y / (FOV/2)
#		def = (strech .+ contract .+ rotate) .* sin.(ωHR*t)
#	end
#	uy(x,y,t) = begin
#		strech = 0 #α * v * (x.^2 .+ y.^2) / (FOV/2)^2 .* sign.(y)
#		contract = - β * v * y / (FOV/2)
#		rotate = γ * v * x / (FOV/2)
#		def = (strech .+ contract .+ rotate) .* sin.(ωHR*t)
#	end
#	# Water spins
#	R = 9/10*FOV/2
#	r = 6/11*FOV/2
#	ring = ⚪(R) .- ⚪(r)
#	ρ = ⚪(r) .+ 0.9*ring #proton density
#	# Diffusion tensor model
#	D = 2e-9 #Diffusion of free water m2/s
#	D1, D2 = D, D/20
#	Dλ1 = D1*⚪(R) #Diffusion map
#	Dλ2 = D1*⚪(r) .+ D2*ring #Diffusion map
#	Dθ =  atan.(x,-y) .* ring #Diffusion map
#	T1 = (1400*⚪(r) .+ 1026*ring)*1e-3 #Myocardial T1
#	T2 = ( 308*⚪(r) .+ 42*ring  )*1e-3 #T2 map [s]
#	# Generating Phantoms
#	heart = Phantom(name="LeftVentricle", x=x, y=y, z=0*x, ρ=0*x, T1=T1, T2=T2, Dλ1=Dλ1, Dλ2=Dλ2, Dθ=Dθ, ux=ux, uy=uy)
#	# Fat spins
#	ring2 = ⚪(FOV/2) .- ⚪(R) #outside fat layer
#	ρ_fat = .5*ρ.*ring2
#	Δw_fat = 2π*220*ring2 #fat should be dependant on B0
#	T1_fat = 800*ring2*1e-3
#	T2_fat = 120*ring2*1e-3 #T2 map [s]
#	fat = Phantom(name="fat", x=x, y=y, z=0*x, ρ=ρ_fat, T1=T1_fat, T2=T2_fat, Δw=Δw_fat)
#	#Resulting phantom
#	obj = fat_bool ? heart + fat : heart #concatenating spins
#end

T, N = 1, 4
seq = RF(1, 1)
seq += Sequence([Grad(1, 1)])
seq += ADC(N, 1)
simParams = KomaMRICore.default_sim_params()
simParams["Δt"], simParams["Δt_rf"] = T/N, T/N
seqd = KomaMRICore.discretize(seq; simParams)
i = Int(floor(length(seqd) / 3))
is_RF_on(seq[1]) == is_RF_on(seqd[1*i:1*i+1]) && is_GR_on(seq[1]) == is_GR_on(seqd[1*i:1*i+1]) && is_ADC_on(seq[1]) == is_ADC_on(seqd[1*i:1*i+1])
is_RF_on(seq[2]) == is_RF_on(seqd[2*i:2*i+1]) && is_GR_on(seq[2]) == is_GR_on(seqd[2*i:2*i+1]) && is_ADC_on(seq[2]) == is_ADC_on(seqd[2*i:2*i+1])
is_RF_on(seq[3]) == is_RF_on(seqd[3*i:3*i+1]) && is_GR_on(seq[3]) == is_GR_on(seqd[3*i:3*i+1]) && is_ADC_on(seq[3]) == is_ADC_on(seqd[3*i:3*i+1])

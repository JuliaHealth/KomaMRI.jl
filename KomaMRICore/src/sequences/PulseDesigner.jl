"""
    PulseDesigner

A module to define different pulse sequences.
"""
module PulseDesigner
using ..KomaMRICore

"""
    ex = RF_hard(B1, T, sys::Scanner; G=[0,0,0], Δf=0)

Definition of the RF hard sequence.

# Arguments
- `B1`: (`Float64`, `[T]`) amplitude of the RF pulse
- `T`: (`Float64`, `[s]`) duration of the RF pulse
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `G`: (`Vector{Float64}`, `=[0, 0, 0]`, `[T]`) gradient amplitudes for x, y, z
- `Δf`: (`Float64`, `=0`, `[Hz]`) frequency offset of the RF pulse

# Returns
- `ex`: (`::Sequence`) excitation Sequence struct

# Examples
```julia-repl
julia> sys = Scanner();

julia> durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse

julia> ex = PulseDesigner.RF_hard(sys.B1, durRF, sys)
Sequence[ τ = 0.587 ms | blocks: 1 | ADC: 0 | GR: 0 | RF: 1 | DEF: 0 ]

julia> plot_seq(ex)
```
"""
RF_hard(B1, T, sys::Scanner; G=[0,0,0], Δf=0) = begin
	ζ = sum(G) / sys.Smax
    gr = reshape([Grad(G[1],T,ζ); Grad(G[2],T,ζ); Grad(G[3],T,ζ)],:,1)
    rf = reshape([RF(B1,T,Δf,ζ)],:,1)
	EX = Sequence(gr, rf)
	EX
end

"""
	seq = spiral_base(FOV::Float64, Nr::Int, sys::Scanner)

Definition of the radial base sequence.
# Arguments
- `FOV`: (`::Float64`, `[m]`) field of view
- `N`: (`::Int`) number of pixels along the radious
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `ex`: (`::Sequence`) RF struct

# References
- MATT A. BERNSTEIN, KEVIN F. KING, XIAOHONG JOE ZHOU, CHAPTER 2 - RADIOFREQUENCY PULSE SHAPES, Handbook of MRI Pulse Sequences, 2004, Pages 35-66, https://doi.org/10.1016/B978-012092861-3/50006-6.
"""
RF_sinc(B1, T, sys::Scanner; G=[0,0,0], Δf=0, a=0.46, TBP=4) = begin
	t0 = T / TBP
	ζ = maximum(abs.(G)) / sys.Smax
	sinc_pulse(t) = B1 * sinc(t/t0) .* ( (1-a) + a*cos((2π*t)/(TBP*t0)) )
	EX = Sequence([	  Grad(G[1],T,ζ) Grad(-G[1],(T-ζ)/2,ζ);	 #Gx
					  Grad(G[2],T,ζ) Grad(-G[2],(T-ζ)/2,ζ);    #Gy
					  Grad(G[3],T,ζ) Grad(-G[3],(T-ζ)/2,ζ)], #Gz
					 [RF(t->sinc_pulse(t - T/2),T; delay=ζ, Δf) RF(0,0)]	 #RF
					)
	EX
end

##################
## Gradient SEQ ##
##################
# # Without moment nulling
# DIF_base(G0,Δ,δ;verbose=false) = begin
# 	DIF = Sequence([Grad(G0,δ) Delay(Δ-δ) -Grad(G0,δ); #Gx
# 	         		Grad(0, δ) Delay(Δ-δ)  Grad(0 ,δ)]) #Gy
# 	if verbose
# 	b_value = get_bvalue(DIF)[1]
# 	GDmax = get_max_grad(DIF)
# 	println("Sequence with diffusion gradient: "*string(round(GDmax*1e3,digits=2))*" mT/m") # s/mm2
# 	println("Sequence with b-value: "*string(round(b_value*1e-6,digits=2))*" s/mm2") # s/mm2
# 	end
# 	DIF
# end
# # With moment nulling: M0,M1,M2
# DIF_null(G0,Δ,δ,verbose=false) = begin
# 	if abs(Δ-4δ/3)<1e-6 #Necessary condition
# 	#Stoecker
# 	DIF = Sequence([Grad(G0,δ/3) -Grad(G0,2δ/3) Delay(Δ-δ) Grad(G0,2δ/3) -Grad(G0,δ/3); #Gx
#              		Delay(δ/3) Delay(2δ/3) Delay(Δ-δ) Delay(2δ/3) Delay(δ/3)]) #Gy
# 	else
# 	#Castillo
# 	DIF = Sequence([ Grad(G0,δ/4) -Grad(G0,δ/2) Grad(G0,δ/4) Delay(Δ-δ)
# 					-Grad(G0,δ/4)  Grad(G0,δ/2)-Grad(G0,δ/4); #Gx
#              		 Delay(δ/4) Delay(δ/2) Delay(δ/4) Delay(Δ-δ)
# 					 Delay(δ/4) Delay(δ/2) Delay(δ/4)])
# 	end

# 	if verbose
# 	b_value = get_bvalue(DIF)[1]
# 	GDmax = get_max_grad(DIF)
# 	println("Sequence with diffusion gradient: "*string(round(GDmax*1e3,digits=2))*" mT/m") # s/mm2
# 	println("Sequence with b-value: "*string(round(b_value*1e-6,digits=2))*" s/mm2") # s/mm2
# 	end
# 	DIF
# end
# EPI
"""
    epi = EPI(FOV::Float64, N::Int, sys::Scanner)

Definition of the EPI sequence.

# Arguments
- `FOV`: (`::Float64`, `[m]`) field of view
- `N`: (`::Int`) number of pixels in the x and y axis
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `epi`: (`::Sequence`) epi Sequence struct

# Examples
```julia-repl
julia> sys, FOV, N = Scanner(), 23e-2, 101

julia> epi = PulseDesigner.EPI(FOV, N, sys)
Sequence[ τ = 62.259 ms | blocks: 203 | ADC: 101 | GR: 205 | RF: 0 | DEF: 4 ]

julia> plot_seq(epi)
```
"""
EPI(FOV::Float64, N::Int, sys::Scanner) = begin
    #TODO: consider when N is even
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	Ga ≥ sys.Gmax ? error("Ga=$(Ga*1e3) mT/m exceeds Gmax=$(Gmax*1e3) mT/m, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	ϵ1 = Δτ/(Δτ+ζ)
	#EPI base
	EPI = Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2],  #Gx
	 	[mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2])) #Gy
	EPI.ADC = [mod(i,2)==1 ? ADC(0,Δτ,ζ) : ADC(N,Ta,ζ) for i=0:2*Ny-2]
	# Relevant parameters
	Δfx_pix = 1/Ta
	Δt_phase = (Ta+2ζ)*Ny + (Δτ+2ζ)*Ny
	Δfx_pix_phase = 1/Δt_phase
	# println("## EPI parameters ##")
	# println("Δx = $(round(Δx*1e3,digits=2)) mm")
	# println("Pixel Δf in freq. direction $(round(Δfx_pix,digits=2)) Hz")
	# println("Pixel Δf in phase direction $(round(Δfx_pix_phase,digits=2)) Hz")
	#Pre-wind and wind gradients
	ϵ2 = Ta/(Ta+ζ)
    PHASE =   Sequence(reshape(1/2*[Grad(      -Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) #This needs to be calculated differently
	DEPHASE = Sequence(reshape(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) #for even N
	seq = PHASE+EPI+DEPHASE
	#Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")
	seq
end

"""
    seq = radial_base(FOV::Float64, Nr::Int, sys::Scanner)

Definition of the radial base sequence.

# Arguments
- `FOV`: (`::Float64`, `[m]`) field of view
- `N`: (`::Int`) number of pixels along the radious
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `seq`: (`::Sequence`) radial base Sequence struct
"""
radial_base(FOV::Float64, Nr::Int, sys::Scanner) = begin
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Δx = FOV/(Nr-1)
	Ta = Δt*(Nr-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	Ga ≥ sys.Gmax ? error("Ga=$(Ga*1e3) mT/m exceeds Gmax=$(Gmax*1e3) mT/m, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	#Radial base
	rad = Sequence([Grad(Ga,Ta,ζ)]) #Gx
	rad.ADC = [ADC(Nr,Ta,ζ)]
	# Acq/Recon parameters
	Nspokes = ceil(Int64, π/2 * Nr ) #Nyquist in the radial direction
	Δθ = π / Nspokes
	# println("## Radial parameters ##")
	# println("FOVr = $(round(FOV*1e2,digits=2)) cm; Δr = $(round(Δx*1e3,digits=2)) mm")
	# println("Nspokes = $Nspokes, to satisfy the Nyquist criterion")
    PHASE = Sequence([Grad(-Ga/2, Ta, ζ)])
	seq = PHASE+rad+PHASE
	#Saving parameters
	seq.DEF = Dict("Nx"=>Nr,"Ny"=>Nspokes,"Nz"=>1,"Δθ"=>Δθ,"Name"=>"radial","FOV"=>[FOV, FOV, 0])

	seq
end

"""
	seq = spiral_base(FOV::Float64, Nr::Int, sys::Scanner)

Definition of the radial base sequence.
# Arguments
- `FOV`: (`::Float64`, `[m]`) field of view
- `N`: (`::Int`) number of pixels along the radious
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `seq`: (`::Function`) function that returns a `Sequence` when evaluated

# References
- Glover, G.H. (1999), Simple analytic spiral K-space algorithm. Magn. Reson. Med., 42: 412-415. https://doi.org/10.1002/(SICI)1522-2594(199908)42:2<412::AID-MRM25>3.0.CO;2-U
"""
spiral_base(FOV::Float64, N::Int64, sys::Scanner; S0=sys.Smax*2/3, Nint=8, λ=Nint/FOV, BW=60e3) = begin
	kmax = N / (2*FOV)
	θmax = kmax / λ # From k(t) = λ θ(t) exp(iθ(t))
	Smax = sys.Smax
	β = Smax * γ / λ
	a₂ = (9β/4)^(1/3)
	Λ = Smax / S0
	θ₁(t) = (.5 * β * t^2) / (Λ + β / 2a₂ * t^(4/3))
	Gmax = sys.Gmax
	ts = (Gmax * 3γ/(2λ*a₂^2))^3 # Gmax = 2λ / 3γ a₂² t ^1/3 e^ (i a₂ t^2/3)
	dt = sys.GR_Δt
	if θ₁(ts) < θmax
		#Region 1 - Slew Rate-Limited
		t1 = 0:dt:ts
		θ₁v = θ₁.(t1)
		#Region 2 - Amplitude-Limtied
		ta = ts .+ (λ/(2γ * sys.Gmax)) * ( θmax.^2 .- θ₁v[end]^2); # ta = ts .+ (λ/(2γ * g0)) * ( θmax.^2 .- θs^2);
		t2 = ts:dt:ta
		θ₂v = sqrt.(θ₁v[end]^2 .+ (2γ/λ) * sys.Gmax * (t2 .- ts))
		θ = [θ₁v[1:end-1]; θ₂v]
	else
		ta = ((2π*FOV)/(3Nint))*sqrt(1/(2γ*Smax*(FOV/N)^3));
		t1 = 0:dt:ta;
		θ = θ₁.(t1)
	end
	dθdt = [0; diff(θ; dims=1) ./ dt]
	#Definition of sequence object
	function spiral(i)
		Δθ = 2π/Nint
		G = (λ/γ) * dθdt .* ( 1 .+ 1im*θ ) .* exp.(1im*(θ .+ Δθ * i))
		Gx = Grad(real.(G),ta,0,abs(real(G[end]))/Smax,0)
		Gy = Grad(imag.(G),ta,0,abs(imag(G[end]))/Smax,0)
		Gz = Grad(0,0)
		GR = reshape([Gx; Gy; Gz], 3, 1)
		R = reshape([RF(0,0)], 1, 1)
		Nadc = floor(Int64, ta*BW)
		A = [ADC(Nadc,ta)]
		seq = Sequence(GR,R,A)
		seq.DEF = Dict("Nx"=>N,"Ny"=>N,"Nz"=>1,"Δθ"=>Δθ,"Nint"=>Nint,"Name"=>"spiral","FOV"=>[FOV, FOV, 0], "λ"=>λ)
		seq
	end
	return spiral
end


export EPI, radial_base
end

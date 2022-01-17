## PulseDesigner
module PulseDesigner
using ..Koma
using ..Koma: γ, Scanner, get_designed_kspace, get_bvalue, get_max_grad

###############
## RF Pulses ##
###############
RF_hard(B1, T, sys::Scanner; G=[0,0,0]) = begin
	ζ = sum(G) / sys.Smax
	EX = Sequence([	Grad(G[1],T,ζ);	 #Gx
					Grad(G[2],T,ζ);  #Gy
					Grad(G[3],T,ζ)], #Gz
					[RF(complex(B1),T,0,ζ)]	 #RF
					)
	EX
end
##################
## Gradient SEQ ##
##################
# Without moment nulling
DIF_base(G0,Δ,δ;verbose=false) = begin
	DIF = Sequence([Grad(G0,δ) Delay(Δ-δ) -Grad(G0,δ); #Gx
	         		Grad(0, δ) Delay(Δ-δ)  Grad(0 ,δ)]) #Gy
	if verbose
	b_value = get_bvalue(DIF)[1]
	GDmax = get_max_grad(DIF)
	println("Sequence with diffusion gradient: "*string(round(GDmax*1e3,digits=2))*" mT/m") # s/mm2
	println("Sequence with b-value: "*string(round(b_value*1e-6,digits=2))*" s/mm2") # s/mm2
	end
	DIF
end
# With moment nulling: M0,M1,M2
DIF_null(G0,Δ,δ,verbose=false) = begin
	if abs(Δ-4δ/3)<1e-6 #Necessary condition
	#Stoecker
	DIF = Sequence([Grad(G0,δ/3) -Grad(G0,2δ/3) Delay(Δ-δ) Grad(G0,2δ/3) -Grad(G0,δ/3); #Gx
             		Delay(δ/3) Delay(2δ/3) Delay(Δ-δ) Delay(2δ/3) Delay(δ/3)]) #Gy
	else
	#Castillo
	DIF = Sequence([ Grad(G0,δ/4) -Grad(G0,δ/2) Grad(G0,δ/4) Delay(Δ-δ)
					-Grad(G0,δ/4)  Grad(G0,δ/2)-Grad(G0,δ/4); #Gx
             		 Delay(δ/4) Delay(δ/2) Delay(δ/4) Delay(Δ-δ)
					 Delay(δ/4) Delay(δ/2) Delay(δ/4)])
	end

	if verbose
	b_value = get_bvalue(DIF)[1]
	GDmax = get_max_grad(DIF)
	println("Sequence with diffusion gradient: "*string(round(GDmax*1e3,digits=2))*" mT/m") # s/mm2
	println("Sequence with b-value: "*string(round(b_value*1e-6,digits=2))*" s/mm2") # s/mm2
	end
	DIF
end
# EPI
EPI_base(FOV::Float64, N::Int, sys::Scanner) = begin
    #TODO: consider when N is odd
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	Ga ≥ sys.Gmax ? error("Error: Ga exceeds Gmax, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	ϵ1 = Δτ/(Δτ+ζ)
	EPI = Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2],#Gx
	 	[mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2])) #Gy
	EPI.ADC = [mod(i,2)==1 ? ADC(0,Δτ,ζ) : ADC(N,Ta,ζ) for i=0:2*Ny-2]
	EPI.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"EPI"=>true)
	# Recon parameters
	try
		Wx = maximum(abs.(get_designed_kspace(EPI)[:,1]))  #no off-set yet so max(k)=Wx
		Wy = maximum(abs.(get_designed_kspace(EPI)[:,2]))
		Δx_pix, Δy_pix = 1/Wx, 1/Wy
		Δfx_pix = γ*Ga*Δx_pix
		Δt_phase = Ta*(Nx/(Nx-1)+(Ny-1)^-1) #Time between two echoes
		Δfx_pix_phase = 1/(Ny*Δt_phase)
		println("## EPI parameters ##")
		println("Δx = "*string(Δx*1e3)*" mm")
		println("Pixel Δf in freq. direction "*string(round(Δfx_pix,digits=2))*" Hz")
		println("Pixel Δf in phase direction "*string(round(Δfx_pix_phase,digits=2))*" Hz")
	catch e
	end
	ϵ2 = Ta/(Ta+ζ)
    PHASE =   Sequence(1/2*[Grad(-Ga, Ta, ζ)      ; ϵ2*Grad(-Ga, Ta, ζ)]) #This needs to be calculated differently
	DEPHASE = Sequence(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)]) #for even N
	seq = PHASE+EPI+DEPHASE
	seq, Δx, Ga, Ta
end

# Radial
radial_base(FOV::Float64, Nr::Int, sys::Scanner) = begin
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Nx = Nr #Number of samples in each spoke
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
	Ga = 1/(γ*Δt*FOV)
	Ga ≥ Gmax ? error("Error: Ga exceeds Gmax, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	rad = Sequence([Grad(Ga,Ta)]) #Gx
	rad.ADC = [ADC(Nr,Ta)]
	# Acq/Recon parameters
	Nspokes = ceil(Int64, π/2 * Nr ) #Nyquist in the radial direction
	Δθ = π / Nspokes
	println("## Radial parameters ##")
	println("FOVr = "*string(round(FOV*1e2,digits=2))*" cm; Δr = "*string(round(Δx*1e3,digits=2))*" mm")
	println("Nspokes = "*string(Nspokes)*", to satisfy the Nyquist criterion")
    PHASE = Sequence([Grad(-Ga, Ta/2)])
	radial = PHASE+rad+PHASE
	radial, Nspokes, Δθ, Ta
end

export DIF_base, DIF_null, EPI_base
end 
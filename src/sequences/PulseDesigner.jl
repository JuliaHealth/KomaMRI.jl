## PulseDesigner
module PulseDesigner
using ..KomaMRI
using ..KomaMRI: γ, Scanner, get_bvalue, get_max_grad

###############
## RF Pulses ##
###############
RF_hard(B1, T, sys::Scanner; G=[0,0,0], Δf=0) = begin
	ζ = sum(G) / sys.Smax
	EX = Sequence([	Grad(G[1],T,ζ);	 #Gx
					Grad(G[2],T,ζ);  #Gy
					Grad(G[3],T,ζ);;], #Gz
					 [RF(B1,T,Δf,ζ);;]	 #RF
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
    PHASE =   Sequence(1/2*[Grad(-Ga, Ta, ζ)      ; ϵ2*Grad(-Ga, Ta, ζ);;]) #This needs to be calculated differently
	DEPHASE = Sequence(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ);;]) #for even N
	seq = PHASE+EPI+DEPHASE
	#Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")
	seq
end

# Radial
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

export EPI, radial_base
end 
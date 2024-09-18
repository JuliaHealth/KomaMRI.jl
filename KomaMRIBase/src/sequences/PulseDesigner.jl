"""
    PulseDesigner

A module to define different pulse sequences.
"""
module PulseDesigner
using ..KomaMRIBase

"""
    seq = RF_hard(B1, T, sys; G=[0, 0, 0], Δf=0)

Returns a sequence with a RF excitation pulse.

# Arguments
- `B1`: (`::Number`, `[T]`) RF pulse amplitude
- `T`: (`::Real`, `[s]`) RF pulse duration
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `G`: (`::Vector{Real}`, `=[0, 0, 0]`, `[T/m]`) gradient amplitudes for x, y, z
- `Δf`: (`::Real`, `=0`, `[Hz]`) RF pulse carrier frequency displacement

# Returns
- `seq`: (`::Sequence`) Sequence struct with a RF pulse

# Examples
```julia-repl
julia> sys = Scanner(); durRF = π / 2 / (2π * γ * sys.B1);

julia> seq = PulseDesigner.RF_hard(sys.B1, durRF, sys);

julia> plot_seq(seq)
```
"""
function RF_hard(B1, T, sys::Scanner; G=[0, 0, 0], Δf=0)
	ζ = sum(G) / sys.Smax
    gr = [Grad(G[1], T, ζ); Grad(G[2], T, ζ); Grad(G[3], T ,ζ) ;;]
    rf = [RF(B1, T, Δf, ζ) ;;]
	return Sequence(gr, rf)
end

"""
	seq = RF_sinc(B1, T, sys; G=[0, 0, 0], Δf=0, a=0.46, TBP=4)

Returns a sequence with a RF sinc waveform.

# References
- Matt A. Bernstein, Kevin F. King, Xiaohong Joe Zhou, Chapter 2 - Radiofrequency Pulse
Shapes, Handbook of MRI Pulse Sequences, 2004, Pages 35-66,
https://doi.org/10.1016/B978-012092861-3/50006-6.

# Arguments
- `B1`: (`::Number`, `[T]`) RF sinc amplitude
- `T`: (`::Real`, `[s]`) RF sinc duration
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `G`: (`::Vector{Real}`, `=[0, 0, 0]`, `[T/m]`) gradient amplitudes for x, y, z
- `Δf`: (`::Real`, `=0`, `[Hz]`) RF pulse carrier frequency displacement
- `a`: (`::Real`, `=0.46`) height appodization window parameter
- `TBP`: (`::Real`, `=4`) width appodization window parameter

# Returns
- `seq`: (`::Sequence`) Sequence struct with a RF pulse

# Examples
```julia-repl
julia> sys = Scanner(); durRF = π / 2 / (2π * γ * sys.B1);

julia> seq = PulseDesigner.RF_sinc(sys.B1, durRF, sys);

julia> plot_seq(seq)
```
"""
function RF_sinc(B1, T, sys::Scanner; G=[0, 0, 0], Δf=0, a=0.46, TBP=4)
	t0 = T / TBP
	ζ = maximum(abs.(G)) / sys.Smax
	sinc_pulse(t) = B1 * sinc(t/t0) .* ((1-a) + a*cos((2π*t)/(TBP*t0)))
    gr1 = [Grad(G[1], T, ζ); Grad(G[2], T, ζ); Grad(G[3], T, ζ)]
    gr2 = [Grad(-G[1], (T-ζ)/2, ζ); Grad(-G[2], (T-ζ)/2, ζ); Grad(-G[3], (T-ζ)/2, ζ)]
    gr = [gr1 gr2]
    rf = [RF(t->sinc_pulse(t - T/2), T; delay=ζ, Δf) RF(0,0)]
	return Sequence(gr, rf)
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
    seq = EPI(FOV::Real, N::Integer, sys::Scanner)

Returns a sequence with EPI gradients.

# Arguments
- `FOV`: (`::Real`, `[m]`) field of view
- `N`: (`::Integer`) number of pixels in the x and y axis
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `seq`: (`::Sequence`) Sequence struct with EPI gradients

# Examples
```julia-repl
julia> sys, FOV, N = Scanner(), 23e-2, 101

julia> seq = PulseDesigner.EPI(FOV, N, sys)

julia> plot_seq(seq)

julia> plot_kspace(seq)
```
"""
function EPI(FOV::Real, N::Integer, sys::Scanner)
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
	GR = zeros(Grad, 3, length(0:2*Ny-2))
	GR.x .= [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2]
	GR.y .= [mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2]
	EPI = Sequence(GR)
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
    PHASE =   Sequence(1/2*[Grad(      -Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ); Grad(0.,0.);;]) #This needs to be calculated differently
	DEPHASE = Sequence(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ); Grad(0.,0.);;]) #for even N
	seq = PHASE+EPI+DEPHASE
	#Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")
	return seq
end

"""
    seq = radial_base(FOV::Real, Nr::Integer, sys::Scanner)

Returns a sequence with radial gradients for a single trajectory.

# Arguments
- `FOV`: (`::Real`, `[m]`) field of view
- `N`: (`::Integer`) number of pixels along the diameter
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `seq`: (`::Sequence`) Sequence struct of a single radial trajectory

# Examples
```julia-repl
julia> sys, FOV, N = Scanner(), 23e-2, 101

julia> seq = PulseDesigner.radial_base(FOV, N, sys)

julia> plot_seq(seq)

julia> plot_kspace(seq)
```
"""
function radial_base(FOV::Real, Nr::Integer, sys::Scanner)
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
	return seq
end

"""
    spiral = spiral_base(FOV, N, sys; S0=sys.Smax*2/3, Nint=8, λ=Nint/FOV, BW=60e3)

Definition of a spiral base sequence.

# References
- Glover, G.H. (1999), Simple analytic spiral K-space algorithm. Magn. Reson. Med.,
42: 412-415. https://doi.org/10.1002/(SICI)1522-2594(199908)42:2<412::AID-MRM25>3.0.CO;2-U

# Arguments
- `FOV`: (`::Real`, `[m]`) field of view
- `N`: (`::Integer`) number of pixels along the radious
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `S0`: (`::Vector{Real}`, `=sys.Smax*2/3`, `[T/m/s]`) slew rate reference
- `Nint`: (`::Integer`, `=8`) number of interleaves
- `λ`: (`::Real`, `=Nint/FOV`, `[1/m]`) kspace spiral parameter
- `BW`: (`::Real`, `=60e3`, `[Hz]`) adquisition parameter

# Returns
- `spiral`: (`::Function`) function that returns a `Sequence` struct when evaluated

# Examples
```julia-repl
julia> sys, FOV, N = Scanner(), 23e-2, 101

julia> spiral = PulseDesigner.spiral_base(FOV, N, sys)

julia> seq = spiral(0)

julia> plot_seq(seq)
```
"""
function spiral_base(
    FOV::Real, N::Integer, sys::Scanner;
    S0=sys.Smax*2/3, Nint=8, λ=Nint/FOV, BW=60e3
)
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
		return seq
	end
	return spiral
end


"""
    seq = EPI_example(; sys=Scanner())

Returns a sequence suitable for acquiring the 2D brain example in the provided examples.

# Keywords
- `sys`: (`::Scanner`) Scanner struct

# Returns
- `seq`: (`::Sequence`) EPI example Sequence struct

# Examples
```julia-repl
julia> seq = PulseDesigner.EPI_example();

julia> plot_seq(seq)
```
"""
function EPI_example(; sys=Scanner())
    B1 = sys.B1;
    durRF = π/2/(2π*γ*B1)
    EX = PulseDesigner.RF_hard(B1, durRF, sys; G=[0,0,0])
    N = 101
    FOV = 23e-2
    EPI = PulseDesigner.EPI(FOV, N, sys)
    TE = 30e-3
    d1 = TE-dur(EPI)/2-dur(EX)
    if d1 > 0 DELAY = Delay(d1) end
    seq = d1 > 0 ? EX + DELAY + EPI : EX + EPI
    seq.DEF["TE"] = round(d1 > 0 ? TE : TE - d1, digits=4)*1e3
    return seq
end

function EPI_dur(duration::Real,FOV::Real,sys::Scanner)
	Δt = sys.ADC_Δt
	FOV_min = 1/(γ*Δt*sys.Gmax)
	@assert FOV ≥ FOV_min "FOV should be at least bigger than $(ceil(FOV_min))"
    let
        epi, epi_new = Sequence(), Sequence()
        for n in 1:2:5001
            epi_new = PulseDesigner.EPI(FOV, n, sys)
            if dur(epi_new) ≤ duration
                epi = epi_new
            else
                d = Delay(duration - dur(epi))
                seq = d + epi
                return seq
            end
        end
    end
end
function EPI_pablo(FOV::Real, N::Integer, sys::Scanner; Δt=sys.ADC_Δt)
	Gmax = sys.Gmax
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	if Ga > Gmax
		println("Ga=$(Ga*1e3) mT/m exceeds Gmax=$(Gmax*1e3) mT/m, increasing Δt to Δt_min="*string(round(1e6/(γ*Gmax*FOV),digits=2))*" us...")
		return EPI_pablo(FOV, N, sys; Δt=1/(γ*Gmax*0.99*FOV))
	end
	ϵ1 = Δτ/(Δτ+ζ)
	#EPI base
	epi = Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2],  #Gx
	 	[mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2])) #Gy
	epi.ADC = [mod(i,2)==1 ? ADC(0,Δτ,ζ) : ADC(N,Ta,ζ) for i=0:2*Ny-2]
	# Relevant parameters
	Δfx_pix = 1/Ta
	Δt_phase = (Ta+2ζ)*Ny + (Δτ+2ζ)*Ny
	Δfx_pix_phase = 1/Δt_phase
	#Pre-wind and wind gradients
	ϵ2 = Ta/(Ta+ζ)
    PHASE =   Sequence(reshape(1/2*[Grad(      -Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) #This needs to be calculated differently
	DEPHASE = Sequence(reshape(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) #for even N
	seq = PHASE+epi+DEPHASE
	#Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")
	return seq
end
function EPI_dif(duration::Real,FOV::Real,sys::Scanner)
	#Δt = sys.ADC_Δt
	#FOV_min = 1/(γ*Δt*sys.Gmax)
	#@assert FOV ≥ FOV_min "FOV should be at least bigger than $(FOV_min)"
	epi = PulseDesigner.EPI_pablo(FOV,128,sys)
	d = Delay(duration - dur(epi))
	seq = d+epi
	return seq
end

"""
    seq = diffusion_SE(b_value::Int, b_vector::Vector, δ::Real, Δ::Real, TE::Real; sys::Scanner = Scanner(), ξ::Real = sys.Gmax/sys.Smax)

Returns a Spin-Echo sequence with two diffusion gradients of width δ and a gap of Δ since the beginning of the first diffusion pulse.

#Arguments
- `b_value`: (`::Int` [s/mm^2]) b-value of the gradient
- `b_vector`: (`::Vector`) normalized position of the gradient in q-space
- `δ`: (`::Real` [s]) duration of each diffusion gradient pulse
- `Δ`: (`::Real` [s]) gap between gradient pulses
- `TE`: (`::Real` [s]) Echo time of the sequence
- `sys`: (`::Scanner`) Scanner struct

#Keywords
- `ξ`: (`::Real` `=sys.Gmax/sys.Smax`) rise/fall times of the gradients

#Returns
- `seq`: (`::Sequence`) Spin-Echo sequence struct with the desired amplitude and gradients' orientation

#Examples
```julia-repl

julia> b_value = 1000

julia> b_vector = [-0.0348169195723329, 0.0869745652700947, -0.995601932052952]

julia> δ = 0.0106

julia> Δ = 0.0735

julia> TE = 150e-3

julia> seq = diffusion_SE(b_value, b_vector, δ, Δ, TE)

julia> plot_seq(seq)
```
"""
function diffusion_SE(b_value::Int, b_vector::Vector, δ::Real, Δ::Real, TE::Real, sys::Scanner; ξ::Real = sys.Gmax/sys.Smax)
	@assert size(b_vector,1) == 3 "b_vector must be a 3 element vector"
	@assert 1 - sqrt(b_vector'*b_vector) < 1e-11 "b_vector must be normalized"
	seq = Sequence()
	b_value < 1e6 ? b_value *=  1e6 : nothing
	b_max = (γ*2*π)^2*sys.Gmax^2*(δ^2*(Δ - δ/3) + ξ^3/3 - δ*ξ^2/6)
	@assert b_value < b_max "With that parameters the maximum allowed b-value is: $(floor(b_max*1e-6))"
	G = Grad(1,δ,ξ,ξ)
	R = hcat(b_vector,zeros(3,2))
	A = sqrt(b_value/b_max)*sys.Gmax
	durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse
	rf = PulseDesigner.RF_hard(sys.B1, durRF, sys)
	rf_inv = (0 + 2.0im)*rf
	delay_grad = Delay((Δ-dur(G))/2)
	@assert TE/2 > durRF/2 + dur(G)+(Δ-dur(G))/2 "With that parameters, the minimum allowed TE is: $((durRF/2+dur(G)+Δ/2)*2e3) ms"
	delay = Delay(TE/2 - durRF - dur(G) - (Δ-dur(G))/2)
	seq = rf+delay+G+delay_grad+rf_inv+delay_grad+G+delay
	@show dur(rf+delay+G+delay_grad+rf_inv+delay_grad+G)*1e3
	seq = R*A*seq
	return seq
end

"""
    seq = diffusion_GRE(b_value::Int, b_vector::Vector, δ::Real, TE::Real,;sys::Scanner = Scanner(), ξ::Real = sys.Gmax/sys.Smax)

Returns a Gradient-Echo sequence with two diffusion gradients of width δ.

#Arguments
- `b_value`: (`::Int` [s/mm^2]) b-value of the gradient
- `b_vector`: (`::Vector`) normalized position of the gradient in q-space
- `δ`: (`::Real` [s]) duration of each diffusion gradient pulse
- `TE`: (`::Real` [s]) Echo time of the sequence
- `sys`: (`::Scanner`) Scanner struct

#Keywords
- `ξ`: (`::Real` `=sys.Gmax/sys.Smax`) rise/fall times of the gradients

#Returns
- `seq`: (`::Sequence`) Spin-Echo sequence struct with the desired amplitude and gradients' orientation

#Examples
```julia-repl

julia> b_value = 100

julia> b_vector = [-0.0348169195723329, 0.0869745652700947, -0.995601932052952]

julia> δ = 0.0106

julia> TE = 50e-3

julia> seq = diffusion_GRE(b_value, b_vector, δ, Δ, TE)

julia> plot_seq(seq)
```
"""
function diffusion_GRE(b_value::Int, b_vector::Vector, δ::Real, TE::Real, sys::Scanner;ξ::Real = sys.Gmax/sys.Smax)
	@assert size(b_vector,1) == 3 "b_vector must be a 3 element vector"
	@assert 1 - sqrt(b_vector'*b_vector) < 1e-11 "b_vector must be normalized"
	seq = Sequence()
	b_value < 1e6 ? b_value = b_value * 1e6 : nothing
	b_max = 2/3*(γ*2*pi)^2*sys.Gmax^2*(δ^2*(2δ/3 + 2ξ)  + ξ^3/3 - δ*ξ/6)
	@assert b_value < b_max "With that parameters the maximum allowed b-value is: $(floor(b_max*1e-6))"
	G = Grad(1,δ,ξ,ξ)
	R = hcat(b_vector,zeros(3,2))
	A = sqrt(b_value/b_max)*sys.Gmax
	durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse
	rf = PulseDesigner.RF_hard(sys.B1, durRF, sys)
	@assert TE/2 > durRF/4 + dur(G) "With that parameters the maximum allowed TE is: $((durRF/4+dur(G))*2e3) ms "	delay = Delay(TE/2 - durRF/4 - dur(G)) #Arbitrary, maybe it should be greater
	delay = Delay(TE/2 - durRF - dur(G))
	seq = rf+delay+G+(-G)+delay
	seq = R*A*seq
	return seq
end


"""
    seq = diffusion_STEAM(b_value::Int, b_vector::Vector, δ::Real, TM::Real, TE::Real,;sys::Scanner = Scanner(), ξ::Real = sys.Gmax/sys.Smax)
Returns a Spin-Echo sequence with two diffusion gradients of width δ and a gap of TM between two last RF pulses.

#Arguments
- `b_value`: (`::Int` [s/mm^2]) b-value of the gradient
- `b_vector`: (`::Vector`) normalized position of the gradient in q-space
- `δ`: (`::Real` [s]) duration of each diffusion gradient pulse
- `TM`: (`::Real` [s]) gap between two last RF pulses
- `TE`: (`::Real` [s]) Echo time of the sequence
- `sys`: (`::Scanner`) Scanner struct

#Keywords
- `ξ`: (`::Real` `=sys.Gmax/sys.Smax`) rise/fall times of the gradients

#Returns
- `seq`: (`::Sequence`) Spin-Echo sequence struct with the desired amplitude and gradients' orientation

#Examples
```julia-repl

julia> b_value = 1000

julia> b_vector = [-0.0348169195723329, 0.0869745652700947, -0.995601932052952]

julia> δ = 0.0106

julia> TM = 0.0735

julia> TE = 150e-3

julia> seq = diffusion_STEAM(b_value, b_vector, δ, TM, TE)

julia> plot_seq(seq)
```
"""
function diffusion_STEAM(b_value::Int, b_vector::Vector, δ::Real, TM::Real, TE::Real, sys::Scanner; ξ::Real = sys.Gmax/sys.Smax)
	@assert size(b_vector,1) == 3 "b_vector must be a 3 element vector"
	@assert 1 - sqrt(b_vector'*b_vector) < 1e-11 "b_vector must be normalized"
	seq = Sequence()
	G = Grad(1,δ,ξ,ξ)
	b_value < 1e6 ? b_value = b_value * 1e6 : nothing
	b_max = (γ*2*pi)^2*sys.Gmax^2*(δ^2*((dur(G)+TM-ξ) - δ/3) + ξ^3/3 - δ*ξ^2/6)
	@assert b_value < b_max "With that parameters the maximum allowed b-value is: $(floor(b_max*1e-6))"
	
	R = hcat(b_vector,zeros(3,2))
	A = sqrt(b_value/b_max)*sys.Gmax
	durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse
	rf = PulseDesigner.RF_hard(sys.B1, durRF, sys)
	@assert TE/2 > durRF/2 + dur(G) "With that parameters, the minimum allowed TE is: $((durRF/2+dur(G)+Δ/2)*2e3) ms"
	delay = Delay(TE/2 - durRF - dur(G))
	seq = rf+delay+G+rf+Delay(TM)+rf+G+delay
	seq = R*A*seq
	return seq
end

function diffusion_SE_EPI(b_value::Int, b_vector::Vector, δ::Real, Δ::Real, TE::Real, FOV::Real, sys::Scanner; ξ::Real = sys.Gmax/sys.Smax)
	@assert size(b_vector,1) == 3 "b_vector must be a 3 element vector"
	# @assert 1 - sqrt(b_vector'*b_vector) < 1e-6 "b_vector must be normalized"
	seq = Sequence()
	b_value < 1e6 ? b_value *=  1e6 : nothing
	b_max = (γ*2*π)^2*sys.Gmax^2*(δ^2*(Δ - δ/3) + ξ^3/3 - δ*ξ^2/6)
	@assert b_value < b_max "With that parameters the maximum allowed b-value is: $(floor(b_max*1e-6))"
	G = Grad(1,δ,ξ,ξ)
	R = hcat(b_vector,zeros(3,2))
	A = sqrt(b_value/b_max)*sys.Gmax
	durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse
	rf = PulseDesigner.RF_hard(sys.B1, durRF, sys)
	rf_inv = (0 + 2.0im)*rf
	delay_grad = Delay((Δ-dur(G))/2)
	@assert TE/2 > durRF/2 + dur(G)+(Δ-dur(G))/2 "With that parameters, the minimum allowed TE is: $((durRF/2+dur(G)+Δ/2)*2e3) ms"
	delay = Delay(TE/2 - durRF - dur(G) - (Δ-dur(G))/2)
	@show epi , t= EPI_dif(2*delay.T,50e-2,sys),2*delay.T
	seq = rf+delay+G+delay_grad+rf_inv+delay_grad+G
	@show dur(rf+delay+G+delay_grad+rf_inv+delay_grad+G)*1e3
	seq = R*A*seq
	return seq+EPI_dif(2*delay.T,FOV,sys)
end

function diffusion_STEAM_EPI(b_value::Int, b_vector::Vector, δ::Real, TM::Real, TE::Real, FOV::Real, sys::Scanner; ξ::Real = sys.Gmax/sys.Smax)
	@assert size(b_vector,1) == 3 "b_vector must be a 3 element vector"
	@assert 1 - sqrt(b_vector'*b_vector) < 1e-5 "b_vector must be normalized"
	seq = Sequence()
	G = Grad(1,δ,ξ,ξ)
	b_value < 1e6 ? b_value = b_value * 1e6 : nothing
	b_max = (γ*2*pi)^2*sys.Gmax^2*(δ^2*((dur(G)+TM-ξ) - δ/3) + ξ^3/3 - δ*ξ^2/6)
	@assert b_value < b_max "With that parameters the maximum allowed b-value is: $(floor(b_max*1e-6))"
	
	R = hcat(b_vector,zeros(3,2))
	A = sqrt(b_value/b_max)*sys.Gmax
	durRF = π/2/(2π*γ*sys.B1); #90-degree hard excitation pulse
	rf = PulseDesigner.RF_hard(sys.B1, durRF, sys)
	@assert TE/2 > durRF/2 + dur(G) "With that parameters, the minimum allowed TE is: $((durRF/2+dur(G)+Δ/2)*2e3) ms"
	delay = Delay(TE/2 - durRF - dur(G))
	seq = rf+delay+G+rf+Delay(TM)+rf+G
	seq = R*A*seq
	return seq+EPI_dif(2*delay.T,FOV,sys)
end

export EPI, radial_base, EPI_example end
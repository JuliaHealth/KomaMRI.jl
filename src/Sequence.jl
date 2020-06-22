#####################
## Sequence OBJECT ##
#####################
# @everywhere begin
"""
    Sequence(GR)

Sequence object.

# Argument
 - `GR::Array{Grad,2}` := Gradient array, fist axis is the dimension (x,y,z) and the second one is time.

# Examples
```julia-repl
julia> A = Sequence([Grad(1,1) for i=1:2,j=1:2])
Sequence(Grad[Grad(1, 1) Grad(1, 1); Grad(1, 1) Grad(1, 1)])

julia> B = Sequence([Grad(2,2) for i=1:2,j=1:2])
Sequence(Grad[Grad(2, 2) Grad(2, 2); Grad(2, 2) Grad(2, 2)])

julia> C = A + B
Sequence(Grad[Grad(1, 1) Grad(1, 1) Grad(2, 2) Grad(2, 2); Grad(1, 1) Grad(1, 1) Grad(2, 2) Grad(2, 2)])
```
"""
struct Sequence
	GR::Array{Grad,2}	#Sequence in (X, Y and Z) and Time
end
# end
Sequence(GR::Array{Grad,1}) = Sequence(reshape(GR,(length(GR),1)))
Sequence() = Sequence([Grad(0,0); Grad(0,0)])
#Sequence operations
+(x::Sequence, y::Sequence) = Sequence([x.GR y.GR])
-(x::Sequence, y::Sequence) = Sequence([x.GR -y.GR])
-(x::Sequence) = Sequence(-x.GR)
*(x::Sequence, α::Real) = Sequence(α*x.GR)
*(α::Real, x::Sequence) = Sequence(α*x.GR)
*(x::Sequence,A::Matrix) = Sequence(x.GR*A)
/(x::Sequence, α::Real) = Sequence(x.GR/α)
#Sequence object functions
is_DAC_on(x::Sequence) = any([g.DAC==1 for g in x.GR])
DAC_on(x::Sequence) = Sequence(DAC_on.(x.GR))
"Duration `T` [s] of the gradient array."
dur(x::Array{Grad,1}) = sum(x[i].T for i=1:size(x,1))
"Duration `T` [s] of the gradient array."
dur(x::Array{Grad,2}) = maximum(sum([x[i,j].T for i=1:size(x,1),j=1:size(x,2)],dims=2))
"Duration `T` [s] of the sequence."
dur(x::Sequence) = dur(x.GR)
"""
	get_grad(seq,dim,t)

Get gradient array from sequence `seq` evaluated in time points `t` for the dimension `dim`.

# Arguments
 - `seq::Sequence`:= Sequence object.
 - `dim::Integer` := Dimension to retrieve, i.e. (1,2,3) = (Gx,Gy,Gz).
 - `t`            := Time vector to evaluate.

# Examples
```julia-repl
julia> seq = Sequence([Grad(1,1) Grad(-1,1); delay(1) delay(1)])
Sequence(Grad[Grad(1, 1) Grad(-1, 1); Grad(0, 1) Grad(0, 1)])

julia> t = 0:.5:2
0.0:0.5:2.0

julia> get_grad(seq,1,t)
5-element Array{Float64,1}:
  1.0
  1.0
  0.0
 -1.0
 -1.0
```
"""
get_grad(seq::Sequence,dim::Integer,t) = begin
	N = size(seq.GR,2)
	A = [seq.GR[dim,i].A for i=1:N]
	T = [seq.GR[dim,i].T for i=1:N]; T = [0; T]
	⊓(t) = (abs.(t.-1/2) .<= 1/2)*1.;
	sum([A[i]*⊓((t.-sum(T[1:i]))/T[i+1]) for i=1:N])
end
get_DAC_on(seq::Sequence,t::Array{Float64,2}) = begin
	Δt = t[2]-t[1]
	M, N = size(seq.GR)
	DAC = sum([seq.GR[i,j].DAC for i=1:M, j=1:N],dims=1)
	T = [seq.GR[1,i].T for i=1:N]; T = [0; T]
	⊓(t) = (abs.(t.-1/2.) .<= 1/2)*1.
	DAC_bool = sum([DAC[i]*⊓((t.-sum(T[1:i]))/(T[i+1]+Δt)) for i=1:N])
	DAC_bool = (DAC_bool.>0)
	circshift(convert.(Bool,DAC_bool[:]),-1)
end
"Calculates the `b`-value, in the PGSE sequnce is b = (2πγGδ)²(Δ-δ/3) and the normalized diffusion vector `n` from a *diffusion sequence*."
get_bvalue(DIF::Sequence) = begin
	M, N = size(DIF.GR)
	G = [DIF.GR[i,j].A for i=1:M,j=1:N] #Strength of pulse
	δ = [DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(DIF) #End of sequence
	# q-vector
	qτ = γ*sum([G[i,j]*δ[i] for i=1:M,j=1:N],dims=2) #Final q-value
	norm(qτ) ≥ 1e-10 ? error("Error: this is not a diffusion sequence, M0{G(t)}=∫G(t)dt!=0.") : 0
	# B-value
	b = zeros(M,N,N)
	for k=1:M,i=1:N,j=1:N
		ij = max(i,j)
		α = i==j ? 2/3 : 1/2
		b[k,i,j] = G[k,i]*G[k,j]*δ[i]*δ[j]*(τ.-T[ij]-α*δ[ij])
	end
	b_value = (2π*γ)^2*sum(b) #Trace of B tensor
	G = [DIF.GR[i,1].A for i=1:M]; A = norm(G)
	b_value, (A!=0) ? G/A : G
end
"Calculates the `b` value for PGSE and a moment nulled sequence."
get_bvalue(G,Δ,δ,null::Bool) = (null ? 1/9 : 1)*(2π*γ*G*δ)^2*(Δ-δ/3)
"Calculates `q` = γδG diffusion vector from a *diffusion sequence*."
get_qvector(DIF::Sequence,type::String) = begin
	M, N = size(DIF.GR)
	G = [DIF.GR[i,j].A for i=1:M,j=1:N] #Strength of pulse
	δ = [DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(DIF) #End of sequence
	# q-vector in time
	qτ = γ*sum([G[i,j]*δ[i] for i=1:M,j=1:N],dims=2) #Final q-value
	norm(qτ) ≥ 1e-10 ? error("Error: this is not a diffusion sequence, M0{G(t)}=∫G(t)dt!=0.") : 0
	# q-trajectory
	int_rect(t,A,δ) = (t.≥0) ? A*(t.-(t.≥δ).*(t.-δ)) : 0
	q(t) = γ*sum([int_rect(t.-T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
 		# Output q-vector and q-trajectory
	if type=="val"
		q(τ/2) #1/m
	elseif type=="traj"
		t = range(0,stop=τ,length=10^3)
		qt = [q(t)[i] for i=1:M, t=t]
		p = plot(t*1e3,qt'*1e-3,label=["qx" "qy" "qz"])
		xlabel!("t [ms]")
		ylabel!("q [1/mm]")
		title!("q(t) trajectory")
		savefig(p,"q-space_traj.pdf")
		t*1e3, qt'*1e-3, p #ms, 1/mm
	end
end
get_M0_M1_M2(SEQ::Sequence,idx=1) = begin
	M, N = size(SEQ.GR)
	G = [SEQ.GR[i,j].A for i=1:M,j=1:N] #Strength of pulse
	δ = [SEQ.GR[idx,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[idx:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(SEQ) #End of sequence
	#M0
	M0i(t,A,δ) = (t.≥0) ? A*(t.-(t.≥δ).*(t.-δ)) : 0
	M0(t) = sum([M0i(t.-T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
	#M1
	M1i(t,A,δ) = (t.≥0) ? A*(t.^2 .-(t.≥δ).*(t.^2 .-δ.^2))./2 : 0
	M1(t) = sum([M1i(t.-T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
	#M1
	M2i(t,A,δ) = (t.≥0) ? A*(t.^3 .-(t.≥δ).*(t.^3 .-δ.^3))./3 : 0
	M2(t) = sum([M2i(t.-T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
	M0, M1, M2
end
get_max_grad(x::Sequence) = begin
	G = [x.GR[i,j].A for i=1:size(x.GR,1),j=1:size(x.GR,2)]
	Gmax = maximum(sqrt.(sum(G.^2,dims=1)))
end
"Outputs designed k-space trajectory from sequence object."
get_designed_kspace(x::Sequence) = begin
	M, N = size(x.GR)
	F = [x.GR[i,j].A*x.GR[i,j].T for j=1:N, i=1:M]
	F = [zeros(1,M); F]
	k = γ*cumsum(F, dims=1)
end
"Outputs actual k-space trajectory for time vector `t` from sequence object."
get_actual_kspace(x::Sequence,t) = begin
	M, N = size(x.GR)
	Nt = length(t)
	Δt = t[2]-t[1]
	G = [get_grad(x,i,t[:]) for i=1:M]
	F = [G[j][i]*Δt for i=1:Nt, j=1:M]
	F = [zeros(1,M); F]
	k = γ*cumsum(F, dims=1)
end
# Without moment nulling
DIF_base(G0,Δ,δ,verbose=false) = begin
	DIF = Sequence([Grad(G0,δ) delay(Δ-δ) -Grad(G0,δ); #Gx
	         		Grad(0, δ) delay(Δ-δ)  Grad(0 ,δ)]) #Gy
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
	DIF = Sequence([Grad(G0,δ/3) -Grad(G0,2δ/3) delay(Δ-δ) Grad(G0,2δ/3) -Grad(G0,δ/3); #Gx
             		delay(δ/3) delay(2δ/3) delay(Δ-δ) delay(2δ/3) delay(δ/3)]) #Gy
	else
	#Castillo
	DIF = Sequence([ Grad(G0,δ/4) -Grad(G0,δ/2) Grad(G0,δ/4) delay(Δ-δ)
					-Grad(G0,δ/4)  Grad(G0,δ/2)-Grad(G0,δ/4); #Gx
             		 delay(δ/4) delay(δ/2) delay(δ/4) delay(Δ-δ)
					 delay(δ/4) delay(δ/2) delay(δ/4)])
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
EPI_base(FOV::Float64,N::Int,Δt::Float64,Gmax::Float64) = begin
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
	Ga = 1/(γ*Δt*FOV)
	Ga ≥ Gmax ? error("Error: Ga exceeds Gmax, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	EPI = DAC_on(Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta) : delay(Ta/(Ny-1)) for i=0:2*Ny-2],#Gx
	 	[mod(i,2)==1 ? Grad(Ga,Ta/(Ny-1)) : delay(Ta) for i=0:2*Ny-2]))) #Gy
	# Recon parameters
	Wx = maximum(abs.(get_designed_kspace(EPI)[:,1]))
	Wy = maximum(abs.(get_designed_kspace(EPI)[:,2]))
	Δx_pix, Δy_pix = 1/Wx, 1/Wy
	Δfx_pix = γ*Ga*Δx_pix
	Δt_phase = Ta*(Nx/(Nx-1)+(Ny-1)^-1) #Time between two echoes
	Δfx_pix_phase = 1/(Ny*Δt_phase)
	println("## EPI parameters ##")
	println("Pixel Δf in freq. direction "*string(round(Δfx_pix,digits=2))*" Hz")
	println("Pixel Δf in phase direction "*string(round(Δfx_pix_phase,digits=2))*" Hz")
    PHASE = Sequence(1/2*[Grad(-Ga, Ta); Grad(-Ga, Ta)])
    DEPHASE = Sequence(1/2*[Grad(Ga, Ta); Grad(-Ga, Ta)])
	PHASE+EPI+DEPHASE, Δx, Ga, Ta
end

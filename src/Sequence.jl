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
mutable struct Sequence
	GR::Array{Grad,2}	#Sequence in (X, Y and Z) and Time
	RF::Array{RF,1}     #RF pulses
	function Sequence(GR) #If no RF is defined, just use a zero amplitude pulse
		new(GR, [RF(0, dur(GR))])
	end
	Sequence(GR,RF) = abs(dur(GR)-dur(RF))>1e-10 ? error("Gradient and RF sequence must have the same duration. 
durGR = $(dur(GR)) s != durRF = $(dur(RF)) s") : new(GR,RF)
end
# end
Sequence(GR::Array{Grad,1}) = Sequence(reshape(GR,(length(GR),1)))
Sequence(GR::Array{Grad,1}, RF::Array{RF,1}) = Sequence(reshape(GR,(length(GR),1)),RF)
Sequence() = Sequence([Grad(0,0); Grad(0,0)])
#Sequence operations
Base.length(x::Sequence) = size(x.GR,2)
Base.iterate(x::Sequence) = (Sequence(x.GR[:,1]), 2)
Base.iterate(x::Sequence, i::Integer) = (i <= length(x)) ? (Sequence(x.GR[:,i]), i+1) : nothing
Base.getindex(x::Sequence, i::UnitRange{Int}) = Sequence(x.GR[:,i])
Base.getindex(x::Sequence, i::Int) = Sequence(x.GR[:,i])
Base.getindex(x::Sequence, i::BitArray{1}) = any(i) ? Sequence(x.GR[:,i]) : nothing
Base.getindex(x::Sequence, i::Array{Bool,1}) = any(i) ? Sequence(x.GR[:,i]) : nothing
Base.lastindex(x::Sequence) = size(x.GR,2)
Base.copy(x::Sequence) where Sequence = Sequence([deepcopy(getfield(x, k)) for k ∈ fieldnames(Sequence)]...)

+(x::Sequence, y::Sequence) = Sequence([x.GR y.GR], [x.RF; y.RF])
-(x::Sequence, y::Sequence) = Sequence([x.GR -y.GR], [x.RF; y.RF])
-(x::Sequence) = Sequence(-x.GR, x.RF)
*(x::Sequence, α::Real) = Sequence(α*x.GR, x.RF)
*(α::Real, x::Sequence) = Sequence(α*x.GR, x.RF)
*(x::Sequence, A::Matrix) = Sequence(x.GR*A, x.RF)
/(x::Sequence, α::Real) = Sequence(x.GR/α, x.RF)
#Sequence object functions
is_DAC_on(x::Sequence) = any([g.DAC for g in x.GR])
"Tells if the sequence has elements with DAC on during t."
is_DAC_on(x::Sequence, t::Array{Float64,2}) = begin
	N = size(x.GR,2)
	T = [i==0 ? 0 : x.GR[1,i].T for i=0:N]
	ts = [sum(T[1:i])  for i=1:N+1]
	activeGrads = [any(ts[i] .<= t .< ts[i+1]) for i=1:N]
	is_DAC_on(x[activeGrads]) #Get elements of sequence active during time period t
end
"Turns on data acquisition and forces RF pulse to have zero amplitud."
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
get_grads(seq::Sequence,t) = begin
	M, N = size(seq.GR);
	A = [seq.GR[i,j].A for i=1:M, j=1:N];
	T = [seq.GR[1,j].T for j=1:N]; T = [0; T];
	⊓(t) = (abs.(t.-1/2) .<= 1/2)*1.;
	(sum([A[j,i]*⊓((t.-sum(T[1:i]))/T[i+1]) for i=1:N]) for j=1:M)
end

"""
Linear interpolation of the gradients, only works for equal gradients' durations.
"""
get_grads_linear(seq::Sequence,t) = begin
	M, N = size(seq.GR);
	A = [seq.GR[i,j].A for i=1:M, j=1:N];
	T = [seq.GR[1,j].T for j=1:N]; T = [0; T];
	∧(t) = max.(1 .- abs.(t),0);
	(sum([A[j,i]*∧((t.-sum(T[1:i]))/T[i+1]) for i=1:N]) for j=1:M)
end

get_DAC_on(seq::Sequence,t::Array{Float64,1}) = begin
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
	G = getproperty.(DIF.GR,:A) #[DIF.GR[i,j].A for i=1:M,j=1:N] #Strength of pulse
	δ = getproperty.(DIF.GR,:T)[1,:] #[DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(DIF) #End of sequence
	# B-value
	b = 0
	for k=1:M,i=1:N,j=1:N
		ij = max(i,j)
		α = (i==j) ? 2/3 : 1/2
		b += G[k,i]*G[k,j]*δ[i]*δ[j]*(τ-T[ij]-α*δ[ij])
	end
	b_value = (2π*γ)^2*b # Trace of B tensor
end
"Calculates the `b`-matrix, such as `b`-value = g' B g [s/mm2] with g [T/m]."
get_Bmatrix(DIF::Sequence) = begin
	M, N = size(DIF.GR)
	δ = getproperty.(DIF.GR,:T)[1,:] #[DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(DIF) #End of sequence
	# B-value
	# b = zeros(M,N,N)
	# for k=1,i=1:N,j=1:N
	# 	ij = max(i,j)
	# 	α = (i==j) ? 2/3 : 1/2
	# 	b[k,i,j] = δ[i]*δ[j]*(τ-T[ij]-α*δ[ij])
	# end
	ij = [max(i,j) for i=1:N, j=1:N]
	α = [(i==j) ? 2/3 : 1/2 for i=1:N, j=1:N]
	b = (δ' .* δ) .* (τ .- T[ij] .- α .* δ[ij])
	b_value = (2π*γ)^2*1e-6*b # Trace of B tensor
end

"Calculates the `b` value for PGSE and a moment nulled sequence."
get_bvalue(G,Δ,δ;null::Bool=false) = (null ? 1/9 : 1)*(2π*γ*G*δ)^2*(Δ-δ/3)
"Calculates `q` = γδG diffusion vector from a *diffusion sequence*."
get_qvector(DIF::Sequence;type::String="val") = begin
	M, N = size(DIF.GR)
	G = getproperty.(DIF.GR,:A) #[DIF.GR[i,j].A for i=1:M,j=1:N] #Strength of pulse
	δ = getproperty.(DIF.GR,:T)[1,:] #[DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(DIF) #End of sequence
	qτ = γ*sum([G[i,j]*δ[j] for i=1:M,j=1:N],dims=2) #Final q-value
	norm(qτ) ≥ 1e-10 ? error("This is not a diffusion sequence, M0{G(t)}=∫G(t)dt!=0.") : 0
	# q-trajectory
	int_rect(t,A,δ) = (t.≥0) ? A*(t.-(t.≥δ).*(t.-δ)) : 0
	q(t) = γ*sum([int_rect(t.-T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
 	# Output q-vector and q-trajectory
	if type=="val"
		q(τ/2) #1/m
	elseif type=="traj"
		t = range(0,stop=τ,length=10^3)
		qt = [q(t)[i] for i=1:M, t=t]
		l = PlotlyJS.Layout(;title="q(t) trajectory", yaxis_title="q [1/mm]",
	    xaxis_title="t [ms]",height=300)
		p = [PlotlyJS.scatter() for j=1:M]
		idx = ["qx" "qy" "qz"]
		for j=1:M
			p[j] = PlotlyJS.scatter(x=t*1e3, y=qt[j,:]*1e-3,name=idx[j])
		end
		PlotlyJS.plot(p, l)
	end
end
get_M0_M1_M2(SEQ::Sequence) = begin
	M, N = size(SEQ.GR)
	G = [SEQ.GR[i,j].A for i=1:M,j=1:N] #Strength of pulse
	δ = [SEQ.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T] #Position of pulse
	τ = dur(SEQ) #End of sequence
	#M0
	M0i(t,T,A,δ) = (t.≥T) ? A*(t.-T.-(t.≥δ+T).*(t.-δ.-T)) : 0
	M0(t) = sum([M0i(t,T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
	#M1
	M1i(t,T,A,δ) = (t.≥T) ? A/2*(t.^2 .-T^2 .-(t.≥δ+T).*(t.^2 .-(δ+T).^2))./2 : 0
	M1(t) = sum([M1i(t,T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
	#M1
	M2i(t,T,A,δ) = (t.≥T) ? A/3*(t.^3 .-T^3 .-(t.≥δ+T).*(t.^3 .-(δ+T).^3))./3 : 0
	M2(t) = sum([M2i(t,T[j],G[i,j],δ[j]) for i=1:M,j=1:N],dims=2)
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
DIF_base(G0,Δ,δ;verbose=false) = begin
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
    #TODO: consider when N is odd
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	Ga ≥ Gmax ? error("Error: Ga exceeds Gmax, increase Δt to at least Δt_min="
	*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	EPI = DAC_on(Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta) : delay(Δτ) for i=0:2*Ny-2],#Gx
	 	[mod(i,2)==1 ? Grad(Ga,Δτ) : delay(Ta) for i=0:2*Ny-2]))) #Gy
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

""""Calculates the normalized moments at the end of the sequence. """
function get_Mmatrix(seq::Sequence)
    M, N = size(seq.GR)
    τ = dur(seq)
    δ = getproperty.(seq.GR,:T)[1,:] #[DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T[1:end-1]] #Position of pulse
    M0 = δ/τ
    M1 = δ.*(T .+ δ/2)/τ^2
	M2 = δ.*(T.^2 .+ T.*δ .+ δ.^2/3)/τ^3
	M3 = δ.*(T.^3 .+ 3/2 * T.^2 .*δ .+ T.*δ.^2 .+ δ.^3/4)/τ^4
    [M0'; M1'; M2'; M3']
end

using LinearAlgebra: I, Bidiagonal, norm
"""Slew rate matrix: |SR*g| ≤ Smax."""
function get_SRmatrix(seq::Sequence)
    _, N = size(seq.GR)
    T = getproperty.(seq.GR[1,:],:T)
    Δt = [(T[i]+T[i+1])/2 for i=1:N-1]
    dv = [Δt; T[end]/2]
    ev = Δt;
	SR = Array(Bidiagonal(-1 ./ dv, 1 ./ ev, :U))
	SR = [SR[1,:]'*0 ; SR]; SR[1,1] = 2/T[1] 
	SR
end

Base.show(io::IO,s::Sequence) = begin
	print(io, "Sequence[ τ = $(round(dur(s)*1e3)) ms ; DAC: $(is_DAC_on(s)) ; GR: $(size(s.GR)) ; RF: $(size(s.RF)) ]")
end

## TO SCANNER
"""Duration in [s] => samples, with dwell-time of Δt = 6.4 μs."""
δ2N(δ) = floor(Int64, δ * 156250) + 2

function write_diff_fwf(DIF, start180, end180, Gmax; filename="./qte_vectors_input.txt")
    open(filename, "w") do io
        dt = 6.4e-6; #6.4us dwell-time
        G1 = DIF[1:start180-1]
        G2 = DIF[end180:end]
        t = range(0,dur(G1),length=δ2N(dur(G1)))
        Gx1, _ = get_grads_linear(G1,t)
        Gx2, _ = get_grads_linear(G2,t)
        write(io, "Maxwell2 $(δ2N(dur(G1))) $(δ2N(dur(G2)))\n")
        for i = 1:length(t)
            fx1 = Gx1[i]/Gmax
            fx2 = -Gx2[i]/Gmax
            line = @sprintf "% .4f % .4f % .4f % .4f % .4f % .4f\n" fx1 fx1 fx1 fx2 fx2 fx2
            write(io, line)
        end
    end
end

function read_diff_fwf(filename="./qte_vectors_input.txt")
    f = readlines(filename)
    G = zeros(size(f[2:end],1),6)
    n1, n2 = tryparse.(Int64, split(f[1]," ")[2:end])
    for (i, l) in enumerate(f[2:end])
        g = tryparse.(Float64, split(l," "))
        G[i,:] = g[g .≠ nothing]
    end
    G, n1, n2
end


#####################
## Sequence OBJECT ##
#####################
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
	GR::Array{Grad,2}		  #Sequence in (X, Y and Z) and time
	RF::Array{RF,2}			  #RF pulses in coil and time
	ADC::Array{ADC,1}		  #ADC in time
	DEF::Dict{String,Any} #Dictionary with information relevant to the reconstructor
	function Sequence(GR)	  #If no RF is DEFined, just use a zero amplitude pulse
		M,N = size(GR)
		new([i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i=1:3, j=1:N], 
			reshape(
				[RF(complex(0), GR[1,i].T, GR[1,i].rise) for i = 1:N]
				,1,N),
			[ADC(0, GR[1,i].T, GR[1,i].rise) for i = 1:N],
			Dict()
			) 
	end
	function Sequence(GR,RF)	#If no ADC is DEFined, just use a ADC with 0 samples
		if size(GR,2) .!= size(RF,2)
			error("The number of Gradient, and RF objects must be the same.")
		end
		M,N = size(GR)
		new([i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i=1:3, j=1:N], 
			RF,
			[ADC(0, GR[1,i].T, GR[1,i].rise) for i = 1:N],
			Dict()) 
	end
	Sequence(GR,RF,ADC) = begin
		if size(GR,2) .!= size(RF,2) .!= length(ADC)
			error("The number of Gradient, RF, and ADC objects must be the same.")
		end
		M,N = size(GR)
		if M != 3
			new([i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i=1:3, j=1:N],RF,ADC, Dict())
		else
			new(GR, RF, ADC, Dict())
		end
	end
	Sequence(GR,RF,ADC,DEF) = begin
		if size(GR,2) .!= size(RF,2) .!= length(ADC)
			error("The number of Gradient, RF, and ADC objects must be the same.")
		end
		M,N = size(GR)
		if M != 3
			new([i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i=1:3, j=1:N],RF,ADC, DEF)
		else
			new(GR, RF, ADC, DEF)
		end
	end
end
#TODO: Add trapezoidal grads MACRO
Sequence(GR::Array{Grad,1}) = Sequence(reshape(GR,(length(GR),1)))
Sequence(GR::Array{Grad,1}, RF::Array{RF,1}) = Sequence(
														reshape(GR,(length(GR),1)),
														reshape(RF,(length(RF),1)),
														[ADC(0,GR[1,i].T) for i = 1:size(GR,2)]
														)
Sequence(GR::Array{Grad,1}, RF::Array{RF,1}, D::ADC, def::Dict) = Sequence(
																reshape(GR,(length(GR),1)),
																reshape(RF,(length(RF),1)),
																[D],
																def
																)
Sequence() = Sequence(reshape([Grad(0,0); Grad(0,0)],2,1))
#Sequence operations. TODO: CHANGE RF to :,i whrn including coil sensitivities
Base.length(x::Sequence) = size(x.GR,2)
Base.iterate(x::Sequence) = (Sequence(x.GR[:,1], x.RF[:,1], x.ADC[1], x.DEF), 2)
Base.iterate(x::Sequence, i::Integer) = (i <= length(x)) ? (Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DEF), i+1) : nothing
Base.getindex(x::Sequence, i::UnitRange{Int}) = Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DEF)
Base.getindex(x::Sequence, i::Int) = Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DEF)
Base.getindex(x::Sequence, i::BitArray{1}) = any(i) ? Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DEF) : nothing
Base.getindex(x::Sequence, i::Array{Bool,1}) = any(i) ? Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DEF) : nothing
Base.lastindex(x::Sequence) = size(x.GR,2)
Base.copy(x::Sequence) where Sequence = Sequence([deepcopy(getfield(x, k)) for k ∈ fieldnames(Sequence)]...)
Base.show(io::IO,s::Sequence) = begin
	compact = get(io, :compact, false)
	if !compact
		print(io, "Sequence[ τ = $(round(dur(s)*1e3)) ms | ADC: $(is_ADC_on(s)) | GR: $(size(s.GR)) | RF: $(size(s.RF)) ]")
	else
		print(io, "Sequence[τ = $(round(dur(s)*1e3)) ms]")
	end
end
#Arithmetic operations
+(x::Sequence, y::Sequence) = Sequence([x.GR y.GR],  [x.RF y.RF], [x.ADC; y.ADC], merge(x.DEF, y.DEF))
-(x::Sequence, y::Sequence) = Sequence([x.GR -y.GR], [x.RF y.RF], [x.ADC; y.ADC], merge(x.DEF, y.DEF))
-(x::Sequence) = Sequence(-x.GR, x.RF, x.ADC, x.DEF)
*(x::Sequence, α::Real) = Sequence(α*x.GR, x.RF, x.ADC, x.DEF)
*(α::Real, x::Sequence) = Sequence(α*x.GR, x.RF, x.ADC, x.DEF)
*(x::Sequence, α::ComplexF64) = Sequence(x.GR, α.*x.RF, x.ADC, x.DEF)
*(α::ComplexF64, x::Sequence) = Sequence(x.GR, α.*x.RF, x.ADC, x.DEF)
*(x::Sequence, A::Matrix{Float64}) = Sequence(A*x.GR, x.RF, x.ADC, x.DEF)
*(A::Matrix{Float64}, x::Sequence) = Sequence(A*x.GR, x.RF, x.ADC, x.DEF)
/(x::Sequence, α::Real) = Sequence(x.GR/α, x.RF, x.ADC, x.DEF)
#Sequence object functions
size(x::Sequence) = size(x.GR[1,:])
is_ADC_on(x::Sequence) = any(x.ADC.N .> 0)
"Tells if the sequence has elements with ADC on during t."
is_ADC_on(x::Sequence, t::Union{Array{Float64,1},Array{Float64,2}}) = begin
	N = length(x)
	ΔT = durs(x)
	ts = cumsum([0 ; ΔT], dims=1)
	activeBlock = [any(ts[i] .<= t .< ts[i+1]) for i=1:N]
	is_ADC_on(x[activeBlock]) #Get elements of sequence active during time period t
end
"Tells if the sequence has elements with RF on."
is_RF_on(x::Sequence) = any(abs.(x.RF.A) .> 0)
"Tells if the sequence has elements with RF on during t."
is_RF_on(x::Sequence, t::Vector{Float64}) = begin
	N = length(x)
	ΔT = durs(x)
	ts = cumsum([0 ; ΔT], dims=1)
	activeBlock = [any(ts[i] .<= t .< ts[i+1]) for i=1:N]
	is_RF_on(x[activeBlock]) #Get elements of sequence active during time period t
end
"Duration of sequence's blocks `ΔT::Array{Real,1}` [s]."
durs(x::Sequence) = begin
	ΔT_GR  = x.GR.dur
	ΔT_RF  = x.RF.dur
	ΔT_ADC = x.ADC.dur
	ΔT = maximum([ΔT_GR ΔT_RF ΔT_ADC],dims=2)
	ΔT
end
"Duration of the sequence `T` [s]."
dur(x::Sequence) = sum(durs(x))

#Trapezoidal Waveform
⏢(t,ΔT,ζ1,ζ2,delay) = begin
	aux = (ζ1    .<= t .- delay .< ζ1+ΔT) .* 1.
	if ζ1 != 0 
		aux .+= (0     .<= t .- delay .< ζ1) .* 1/2 #(t .- delay)./ζ1 
	end
	if ζ2 !=0
		aux .+= (ζ1+ΔT .<= t .- delay .< ζ1+ΔT+ζ2) .* 1/2 #(t.-delay.-ζ1.-ΔT)./ζ2
	end
	aux
end
"""
	get_grad(seq,dim,t)

Get gradient array from sequence `seq` evaluated in time points `t` for the dimension `dim`.

# Arguments
 - `seq::Sequence`:= Sequence object.
 - `dim::Integer` := Dimension to retrieve, i.e. (1,2,3) = (Gx,Gy,Gz).
 - `t`            := Time vector to evaluate.

# Examples
```julia-repl
julia> seq = Sequence([Grad(1,1) Grad(-1,1); Delay(1) Delay(1)])
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
	#Amplitude
	A = seq.GR.A
	#Timings
	T = seq.GR.T
	ζ1 = seq.GR.rise
	ζ2 = seq.GR.fall
	delay = seq.GR.delay
	T0 = cumsum([0; durs(seq)], dims=1)
	(sum([A[j,i]*⏢(t.-T0[i],T[i],ζ1[i],ζ2[i],delay[i]) for i=1:length(seq)]) for j=1:3)
end

get_rfs(seq::Sequence,t) = begin
	#Amplitude
	A = seq.RF.A
	#Timings
	T = seq.RF.T
	delay = seq.RF.delay
	T0 = cumsum([0; durs(seq)], dims=1)
	[sum([A[1,i]*⏢(t.-T0[i],T[i],0,0,delay[i]) for i=1:length(seq)])]
end

"""
Linear interpolation of the gradients, only works for gradients with uniform duration.
"""
get_grads_linear(seq::Sequence,t) = begin
	M, N = size(seq.GR);
	A = [seq.GR[i,j].A for i=1:M, j=1:N];
	T = [seq.GR[1,j].T for j=1:N]; T = [0; T];
	∧(t) = max.(1 .- abs.(t),0);
	(sum([A[j,i]*∧((t.-sum(T[1:i]))/T[i+1]) for i=1:N]) for j=1:M)
end

get_ADC_on(seq::Sequence,t::Array{Float64,1}) = begin
	Δt = t[2]-t[1]
	M, N = size(seq.ADC)
	ADC = seq.ADC.N .> 0
	T = [seq.ADC[i].T for i=1:N]; T = [0; T]
	⊓(t) = (0 .<= t .< 1)*1.;
	ADC_bool = sum([ADC[i]*⊓((t.-sum(T[1:i]))/(T[i+1]+Δt)) for i=1:N])
	ADC_bool = (ADC_bool.>0)
	circshift(convert.(Bool,	[:]),-1)
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
	# B-value, slower way of doing it
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
	b_value
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
	x = x[(!).(Koma.is_RF_on.(x))] #Remove RF blocks
	M, N = size(x.GR)
	rise = x.GR.A .* x.GR.rise' / 2.
	flat = x.GR.A .* x.GR.T'
	fall = x.GR.A .* x.GR.fall' / 2.
	F = cat(([rise[:,i] flat[:,i] fall[:,i]] for i=1:N)...; dims=2)
	F = [zeros(1,M); F']
	k = γ*cumsum(F, dims=1)
	k
end
"Outputs actual k-space trajectory for time vector `t` from sequence object."
get_actual_kspace(x::Sequence, simsParams) = begin
	seq = x[(!).(is_RF_on.(x))] #Remove RF blocks
	t, Δt = get_uniform_times(seq, simsParams["Δt"])
	Gx, Gy, Gz = get_grads(seq,t)
	G = cat(Gx,Gy,Gz; dims=2)
	F = G .* Δt 
	F = [zeros(1,3); F]
	k = γ*cumsum(F, dims=1)
	k
end

""""Calculates the normalized moments at the end of the sequence. """
function get_Mmatrix(seq::Sequence;order=0)
    M, N = size(seq.GR)
    τ = dur(seq)
    δ = getproperty.(seq.GR,:T)[1,:] #[DIF.GR[1,j].T for j=1:N]; #Duration of pulse
	T = [sum(δ[1:j]) for j=1:N]; T = [0; T[1:end-1]] #Position of pulse
	if order == 0
		M0 = δ/τ
		M1 = δ.*(T .+ δ/2)/τ^2
		M2 = δ.*(T.^2 .+ T.*δ .+ δ.^2/3)/τ^3
		M3 = δ.*(T.^3 .+ 3/2 * T.^2 .*δ .+ T.*δ.^2 .+ δ.^3/4)/τ^4
	elseif order == 1
		M0 = δ/τ
		M1 = δ.*T/τ^2
		M2 = δ.*(T.^2 .+ δ.^2/6)/τ^3
		M3 = δ.*(T.^3 .+ T.*δ.^2/2)/τ^4
	end
    [M0'; M1'; M2'; M3']
end

using LinearAlgebra: I, Bidiagonal, norm
"""Slew rate matrix: |SR*g| ≤ Smax."""
function get_SRmatrix(seq::Sequence)
    _, N = size(seq.GR)
    T = getproperty.(seq.GR[1,:],:T)
    Δt = [T[i] for i=1:N-1]
    dv = [Δt; T[end]]
    ev = Δt;
	SR = Array(Bidiagonal(-1 ./ dv, 1 ./ ev, :U))
	SR = [SR[1,:]'*0 ; SR]; SR[1,1] = T[1] 
	SR
end

## TO SCANNER (Philips)
"""Duration in [s] => samples, with dwell-time of Δt = 6.4 μs."""
δ2N(δ) = floor(Int64, δ * 156250) + 2

function write_diff_fwf(DIF, idx180, Gmax, bmax; filename="./qte_vectors_input.txt", name="Maxwell2")
    open(filename, "w") do io
        G1 = DIF[1:idx180-1]
        G2 = DIF[idx180 + 1:end]
        t1 = range(0,dur(G1),length=δ2N(dur(G1)))
		t2 = range(0,dur(G2),length=δ2N(dur(G2)))
		Gx2, Gy2, Gz2 = get_grads_linear(G2,t2)
		Gx1, Gy1, Gz1 = get_grads_linear(G1,t1)
        write(io, "$name $(δ2N(dur(G1))) $(δ2N(dur(G2))) $(bmax)\n")
        for i = 1:length(t1)
            fx1, fy1, fz1 = (Gx1[i], Gy1[i], Gz1[i])./Gmax
            fx2, fy2, fz2 = i ≤ length(t2) ? (-Gx2[i], -Gy2[i], -Gz2[i])./Gmax : (0,0,0)
            line = @sprintf "% .4f % .4f % .4f % .4f % .4f % .4f\n" fx1 fy1 fz1 fx2 fy2 fz2
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


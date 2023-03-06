"""
    seq = Sequence()
    seq = Sequence(GR)
    seq = Sequence(GR, RF)
    seq = Sequence(GR, RF, ADC)
    seq = Sequence(GR, RF, ADC, DUR)
    seq = Sequence(GR::Array{Grad,1})
    seq = Sequence(GR::Array{Grad,1}, RF::Array{RF,1})
    seq = Sequence(GR::Array{Grad,1}, RF::Array{RF,1}, A::ADC, DUR, DEF)

The Sequence struct.

# Arguments
- `GR`: (`::Matrix{Grad}`) gradient matrix, rows are for (x,y,z) and columns are for time
- `RF`: (`::Matrix{RF}`) RF matrix, the 1 row is for the coil and columns are for time
- `ADC`: (`::Vector{ADC}`) ADC vector in time
- `DUR`: (`::Vector{Float64}`, `[s]`) duration of each sequence-block, this enables
    delays after RF pulses to satisfy ring-down times
- `DEF`: (`::Dict{String, Any}`) dictionary with relevant information of the sequence.
    The possible keys are [`"AdcRasterTime"`, `"GradientRasterTime"`, `"Name"`, `"Nz"`,
    `"Num_Blocks"`, `"Nx"`, `"Ny"`, `"PulseqVersion"`, `"BlockDurationRaster"`,
    `"FileName"`, `"RadiofrequencyRasterTime"`]

# Returns
- `seq`: (`::Sequence`) Sequence struct
"""
mutable struct Sequence
	GR::Array{Grad,2}		  #Sequence in (X, Y and Z) and time
	RF::Array{RF,2}			  #RF pulses in coil and time
	ADC::Array{ADC,1}		  #ADC in time
	DUR::Vector				  #Duration of each block, this enables delays after RF pulses to satisfy ring-down times
	DEF::Dict{String,Any} 	  #Dictionary with information relevant to the reconstructor
end

#MAIN CONSTRUCTORS
function Sequence(GR)	  #If no RF is defined, just use a zero amplitude pulse
	M,N = size(GR)
	Sequence([i <= M ? GR[i,j] : Grad(0, 0) for i=1:3, j=1:N],
		reshape([RF(0, 0) for i = 1:N],1,:),
		[ADC(0, 0) for i = 1:N],
		GR.dur,
		Dict()
		)
end
function Sequence(GR, RF)	#If no ADC is DEFined, just use a ADC with 0 samples
	@assert size(GR,2) .== size(RF,2) "The number of Gradient, and RF objects must be the same."
	M,N = size(GR)
	Sequence([i <= M ? GR[i,j] : Grad(0, 0) for i=1:3, j=1:N],
		RF,
		[ADC(0, 0) for i = 1:N],
		maximum([GR.dur RF.dur],dims=2)[:],
		Dict()
		)
end
function Sequence(GR, RF, ADC)
	@assert size(GR,2) .== size(RF,2) .== length(ADC) "The number of Gradient, RF, and ADC objects must be the same."
	M,N = size(GR)
	Sequence([i <= M ? GR[i,j] : Grad(0, 0) for i=1:3, j=1:N],
		RF,
		ADC,
		maximum([GR.dur RF.dur ADC.dur],dims=2)[:],
		Dict())
end
function Sequence(GR, RF, ADC, DUR)
	@assert size(GR,2) .== size(RF,2) .== length(ADC) .== length(DUR) "The number of Gradient, RF, ADC, and DUR objects must be the same."
	M,N = size(GR)
	Sequence([i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i=1:3, j=1:N],
		RF,
		ADC,
		maximum([GR.dur RF.dur ADC.dur DUR],dims=2)[:],
		Dict())
end

#OTHER CONSTRUCTORS
Sequence(GR::Array{Grad,1}) = Sequence(reshape(GR,1,:))
Sequence(GR::Array{Grad,1}, RF::Array{RF,1}) = Sequence(
														reshape(GR,:,1),
														reshape(RF,1,:),
														[ADC(0,0) for i = 1:size(GR,2)]
														)
Sequence(GR::Array{Grad,1}, RF::Array{RF,1}, A::ADC, DUR, DEF) = Sequence(
																reshape(GR,:,1),
																reshape(RF,1,:),
																[A],
																[DUR],
																DEF
																)
Sequence() = Sequence([Grad(0,0)])

"""
    str = show(io::IO, s::Sequence)

Displays information about the Sequence struct `s` in the julia REPL.

# Arguments
- `s`: (`::Sequence`) Sequence struct

# Returns
- `str` (`::String`) output string message
"""
Base.show(io::IO, s::Sequence) = begin
	compact = get(io, :compact, false)
	if !compact
		nGRs = sum(is_Gx_on.(s)) + sum(is_Gy_on.(s)) + sum(is_Gz_on.(s))
		print(io, "Sequence[ τ = $(round(dur(s)*1e3;digits=3)) ms | blocks: $(length(s)) | ADC: $(sum(is_ADC_on.(s))) | GR: $nGRs | RF: $(sum(is_RF_on.(s))) | DEF: $(length(s.DEF)) ]")
	else
		print(io, "Sequence[τ = $(round(dur(s)*1e3;digits=3)) ms]")
	end
end

#Sequence operations
Base.length(x::Sequence) = length(x.DUR)
Base.iterate(x::Sequence) = (Sequence(x.GR[:,1], x.RF[:,1], x.ADC[1], x.DUR[1], x.DEF), 2)
Base.iterate(x::Sequence, i::Integer) = (i <= length(x)) ? (Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DUR[i], x.DEF), i+1) : nothing
Base.getindex(x::Sequence, i::UnitRange{Int}) = Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DUR[i], x.DEF)
Base.getindex(x::Sequence, i::Int) = Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DUR[i], x.DEF)
Base.getindex(x::Sequence, i::BitArray{1}) = any(i) ? Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DUR[i], x.DEF) : nothing
Base.getindex(x::Sequence, i::Array{Bool,1}) = any(i) ? Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DUR[i], x.DEF) : nothing
Base.lastindex(x::Sequence) = length(x.DUR)
Base.copy(x::Sequence) where Sequence = Sequence([deepcopy(getfield(x, k)) for k ∈ fieldnames(Sequence)]...)

#Arithmetic operations
+(x::Sequence, y::Sequence) = Sequence([x.GR  y.GR], [x.RF y.RF], [x.ADC; y.ADC], [x.DUR; y.DUR], merge(x.DEF, y.DEF))
-(x::Sequence, y::Sequence) = Sequence([x.GR -y.GR], [x.RF y.RF], [x.ADC; y.ADC], [x.DUR; y.DUR], merge(x.DEF, y.DEF))
-(x::Sequence) = Sequence(-x.GR, x.RF, x.ADC, x.DUR, x.DEF)
*(x::Sequence, α::Real) = Sequence(α*x.GR, x.RF, x.ADC, x.DUR, x.DEF)
*(α::Real, x::Sequence) = Sequence(α*x.GR, x.RF, x.ADC, x.DUR, x.DEF)
*(x::Sequence, α::ComplexF64) = Sequence(x.GR, α.*x.RF, x.ADC, x.DUR, x.DEF)
*(α::ComplexF64, x::Sequence) = Sequence(x.GR, α.*x.RF, x.ADC, x.DUR, x.DEF)
*(x::Sequence, A::Matrix{Float64}) = Sequence(A*x.GR, x.RF, x.ADC, x.DUR, x.DEF) #TODO: change this, Rotation fo waveforms is broken
*(A::Matrix{Float64}, x::Sequence) = Sequence(A*x.GR, x.RF, x.ADC, x.DUR, x.DEF) #TODO: change this, Rotation fo waveforms is broken
/(x::Sequence, α::Real) = Sequence(x.GR/α, x.RF, x.ADC, x.DUR, x.DEF)
#Grad operations
+(s::Sequence, g::Grad) = s + Sequence([g;;])
+(g::Grad, s::Sequence) = Sequence([g;;]) + s
#RF operations
+(s::Sequence, r::RF) = s + Sequence([Grad(0,0);;],[r;;])
+(r::RF, s::Sequence) = Sequence([Grad(0,0);;],[r;;]) + s
#ADC operations
+(s::Sequence, a::ADC) = s + Sequence([Grad(0,0);;],[RF(0,0);;],[a])
+(a::ADC, s::Sequence) = Sequence([Grad(0,0);;],[RF(0,0);;],[a]) + s
#Sequence object functions
size(x::Sequence) = size(x.GR[1,:])

"""
    y = is_ADC_on(x::Sequence)
    y = is_ADC_on(x::Sequence, t::Union{Array{Float64,1}, Array{Float64,2}})

Tells if the sequence `seq` has elements with ADC active, or active during time `t`.

# Arguments
- `x`: (`::Sequence`) sequence struct
- `t`: (`::Union{Array{Float64,1}, Array{Float64,2}}`, `[s]`) time to check

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the ADC in the sequence is active
"""
is_ADC_on(x::Sequence) = any(x.ADC.N .> 0)
is_ADC_on(x::Sequence, t::AbstractVecOrMat) = begin
	N = length(x)
	ΔT = durs(x)
	ts = cumsum([0 ; ΔT[1:end-1]])
	delays = x.ADC.delay
	Ts = 	 x.ADC.dur #delat+T
	t0s = ts .+ delays
	tfs = ts .+ Ts
	# The following just checks the ADC
	# |___∿  |
	#     △
	#     Here
	activeADC = any([is_ADC_on(x[i]) && any(t0s[i] .<= t .< tfs[i]) for i=1:N])
	activeADC
end

"""
    y = is_RF_on(x::Sequence)
    y = is_RF_on(x::Sequence, t::Vector{Float64})

Tells if the sequence `seq` has elements with RF active, or active during time `t`.

# Arguments
- `x`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) time to check

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the RF in the sequence is active
"""
is_RF_on(x::Sequence) = any([sum(abs.(r.A)) for r = x.RF] .> 0)
is_RF_on(x::Sequence, t::AbstractVector) = begin
	N = length(x)
	ΔT = durs(x)
	ts = cumsum([0 ; ΔT[1:end-1]])
	delays = x.RF.delay
	Ts = 	 x.RF.dur #dur = delat+T
	t0s = ts .+ delays
	tfs = ts .+ Ts
	# The following just checks the RF waveform
	# |___∿  |
	#     △
	#     Here
	activeRFs = any([is_RF_on(x[i]) && any(t0s[i] .<= t .<= tfs[i]) for i=1:N])
	activeRFs
end

"""
    y = is_GR_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GR in the sequence is active
"""
is_GR_on(x::Sequence) = any([sum(abs.(r.A)) for r = x.GR] .> 0)

"""
    y = is_Gx_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active in x direction.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GRx in the sequence is active
"""
is_Gx_on(x::Sequence) = any([sum(abs.(r.A)) for r = x.GR.x] .> 0)

"""
    y = is_Gy_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active in y direction.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GRy in the sequence is active
"""
is_Gy_on(x::Sequence) = any([sum(abs.(r.A)) for r = x.GR.y] .> 0)

"""
    y = is_Gz_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active in z direction.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GRz in the sequence is active
"""
is_Gz_on(x::Sequence) = any([sum(abs.(r.A)) for r = x.GR.z] .> 0)

"""
    y = is_Delay(x::Sequence)

Tells if the sequence `seq` is a delay.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y::Bool`: boolean that tells whether or not the sequence is a delay
"""
is_Delay(x::Sequence) = !(is_GR_on(x) || is_RF_on(x) || is_ADC_on(x))

"""
    ΔT = durs(x::Sequence)

Returns the array of durations of sequence's blocks in [s].

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `ΔT`: (`::Vector{Real}`, `[s]`) array of durations of sequence's blocks
"""
durs(x::Sequence) = begin
	# ΔT_GR  = x.GR.dur
	# ΔT_RF  = x.RF.dur
	# ΔT_ADC = x.ADC.dur
	# ΔT_DUR = x.DUR
	# ΔT = maximum([ΔT_GR ΔT_RF ΔT_ADC ΔT_DUR],dims=2)
	x.DUR
end

"""
    T = dur(x::Sequence)

The total duration of the sequence in [s].

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `T`: (`::Real`, `[s]`) total duration of the sequence
"""
dur(x::Sequence) = sum(durs(x))

"""
    y = ⏢(A, t, ΔT, ζ1, ζ2, delay)

Generates a trapezoidal waveform vector.

# Arguments
- `A`: (`::Real`) amplitude
- `t`: (`::Vector{Float64}`, `[s]`) times to evaluate (actually it's a `1-row
    ::Matrix{Float64}`)
- `ΔT`: (`::Real`, `[s]`) time duration of the top-flat
- `ζ1`: (`::Real`, `[s]`) rise time duration
- `ζ2`: (`::Real`, `[s]`) fall time duration
- `delay`: (`::Real`, `[s]`) delay time

# Returns
- `y`: (`::Vector{Float64}`) trapezoidal waveform (actually it's a `1-row
    ::Matrix{Float64}`)
"""
⏢(A, t, ΔT, ζ1, ζ2, delay) = begin
	if sum(abs.(A)) != 0 && ΔT+ζ1+ζ2 != 0 # If no event just ignore calculations
		#Getting amplitudes, only supports uniformly sampled waveforms for now
		if length(A) != 1
			grad_raster = ΔT / length(A)
			idx = ceil.(Int, (t .- delay .- ζ1) ./ grad_raster ) #Time to integer index
			valid = 1 .<= idx .<= length(A)
			idx[(!).(valid)] .= 1
			B = A[idx] .* valid
		else
			B = A
		end
		#Trapezoidal waveform
		aux = (ζ1    .<= t .- delay .< ζ1+ΔT) .* B
		if ζ1 != 0
			aux .+= (0     .<= t .- delay .< ζ1) .* A[1] .* (t .- delay)./ζ1
		end
		if ζ2 !=0
			aux .+= (ζ1+ΔT .<= t .- delay .< ζ1+ΔT+ζ2) .* A[end] .* (1 .- (t.-delay.-ζ1.-ΔT)./ζ2)
		end
	else
		aux = zeros(size(t))
	end
	aux
end

"""
    Gx, Gy, Gz = get_grads(seq, t::Vector)
    Gx, Gy, Gz = get_grads(seq, t::Matrix)

Get the gradient array from sequence `seq` evaluated in time points `t`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}` or `1-row ::Matrix{Float64}`, `[s]`) times to evaluate

# Returns
- `Gx`: (`Vector{Float64}` or `1-row ::Matrix{Float64}`, `[T]`) gradient vector values
    in the x direction
- `Gy`: (`Vector{Float64}` or `1-row ::Matrix{Float64}`, `[T]`) gradient vector values
    in the y direction
- `Gz`: (`Vector{Float64}` or `1-row ::Matrix{Float64}`, `[T]`) gradient vector values
    in the z direction
"""
function get_grads(seq, t::Vector)
    gx = get_theo_Gi(seq, 1)
    gy = get_theo_Gi(seq, 2)
    gz = get_theo_Gi(seq, 3)
    Gx = linear_interpolation(gx..., extrapolation_bc=0)(t)
    Gy = linear_interpolation(gy..., extrapolation_bc=0)(t)
    Gz = linear_interpolation(gz..., extrapolation_bc=0)(t)
    (Gx, Gy, Gz)
end
function get_grads(seq, t::Matrix)
	t_vec = t[:]
    gx = get_theo_Gi(seq, 1)
    gy = get_theo_Gi(seq, 2)
    gz = get_theo_Gi(seq, 3)
    Gx = linear_interpolation(gx..., extrapolation_bc=0)(t_vec)
    Gy = linear_interpolation(gy..., extrapolation_bc=0)(t_vec)
    Gz = linear_interpolation(gz..., extrapolation_bc=0)(t_vec)
    (Gx', Gy', Gz')
end
# hold_interpolation(range::AbstractVector, vs::AbstractVector; extrapolation_bc = Throw()) =
#     extrapolate(interpolate((range, ), vs, Gridded(Constant{Previous}())), 0)
# get_grads(seq::Sequence,t) = begin
# 	#Amplitude
# 	A = seq.GR.A
# 	#Grad Timings
# 	T = seq.GR.T
# 	ζ1 = seq.GR.rise
# 	ζ2 = seq.GR.fall
# 	delay = seq.GR.delay
# 	#Sequence timings
# 	ΔT = durs(seq) #Duration of sequence block
# 	T0 = cumsum([0; ΔT[:]]) #Start time of each block
# 	#Waveforms
# 	(sum([⏢(A[j,i],t.-T0[i],T[j,i],ζ1[j,i],ζ2[j,i],delay[j,i]) for i=1:length(seq)]) for j=1:3)
# end

"""
    B1, Δf_rf  = get_rfs(seq::Sequence, t)

Returns the RF pulses and the delta frequency.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`1-row ::Matrix{Float64}`, `[s]`) time points

# Returns
- `B1`: (`1-row ::Matrix{ComplexF64}`, `[T]`) vector of RF pulses
- `Δf_rf`: (`1-row ::Matrix{Float64}`, `[Hz]`) delta frequency vector
"""
get_rfs(seq::Sequence, t) = begin
	#Amplitude
	A  = seq.RF.A
	Δf = seq.RF.Δf
	#Timings
	T = seq.RF.T
	delay = seq.RF.delay
	T0 = cumsum([0; durs(seq)], dims=1)
	(sum([⏢(A[1,i], t.-T0[i],sum(T[i]),0,0,delay[i]) for i=1:length(seq)]),
	 sum([⏢(Δf[1,i],t.-T0[i],sum(T[i]),0,0,delay[i]) for i=1:length(seq)])
	)
end
"""
    y = get_flip_angles(x::Sequence)

Returns all the flip angles of the RF pulses in the sequence `x`.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Vector{Float64}`, `[deg]`) flip angles
"""
get_flip_angles(x::Sequence) = get_flip_angle.(x.RF)[:]

"""
    y = get_ADC_on(seq::Sequence, t::Array{Float64,1})

Get the ADC struct that are active.

!!! note
    This function is not being used.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) the time vector

# Returns
- `y`: (`::Vector{Bool}`) ADC struct that are active
"""
get_ADC_on(seq::Sequence, t::Array{Float64,1}) = begin
	Δt = t[2]-t[1]
	M, N = size(seq.ADC)
	ADC = seq.ADC.N .> 0
	T = [seq.ADC[i].T for i=1:N]; T = [0; T]
	⊓(t) = (0 .<= t .< 1)*1.;
	ADC_bool = sum([ADC[i]*⊓((t.-sum(T[1:i]))/(T[i+1]+Δt)) for i=1:N])
	ADC_bool = (ADC_bool.>0)
	circshift(convert.(Bool,	[:]),-1)
end

"""
    bvalue = get_bvalue(DIF::Sequence)

Calculates the `b`-value, in the PGSE sequnce is b = (2πγGδ)²(Δ-δ/3) and the normalized
diffusion vector `n` from a *diffusion sequence*.

    bvalue = get_bvalue(G, Δ, δ; null::Bool=false)

Calculates the `b` value for PGSE and a moment nulled sequence.

!!! note
    This function is not being used.

# Arguments
- `DIF`: (`:Sequence`)
- `G`
- `Δ`
- `δ`

# Keywords
- `null`: (`::Bool`, `=false`)

# Returns
- `bvalue`: the b value
"""
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

get_bvalue(G, Δ, δ; null::Bool=false) = (null ? 1/9 : 1)*(2π*γ*G*δ)^2*(Δ-δ/3)

"""
    bmatrix = get_Bmatrix(DIF::Sequence)

Calculates the `b`-matrix, such as `b`-value = g' B g [s/mm2] with g [T/m].

!!! note
    This function is not being used.

# Arguments
- `DIF`: (`::Sequence`) diffusion Sequence struct

# Returns
- `bmatrix`: b matrix
"""
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

"""
    qvector = get_qvector(DIF::Sequence; type::String="val")

Calculates `q` = γδG diffusion vector from a *diffusion sequence*.

!!! note
    This function is not being used.

# Arguments
- `DIF`: (`::Sequence`) diffusion Sequence struct
- `type`: (`::String`, `="val"`) option to get a value or plot

# Returns
- `qvector`: q diffusion vector as a value or a plot
"""
get_qvector(DIF::Sequence; type::String="val") = begin
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

"""
    M0, M1, M2 = get_M0_M1_M2(seq::Sequence)

Get the momentums `M0`, `M1` and `M2` from the sequence `seq`.

!!! note
    This function is not being used.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `M0`: M0 momentum
- `M1`: M1 momentum
- `M2`: M2 momentum
"""
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

"""
    Gmax = get_max_grad(x::Sequence)

The maximum gradient of the sequence `x`.

!!! note
    This function is not being used.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `Gmax`: (`::Real`) maximum value of the gradient
"""
get_max_grad(x::Sequence) = begin
	G = [x.GR[i,j].A for i=1:size(x.GR,1),j=1:size(x.GR,2)]
	Gmax = maximum(sqrt.(sum(G.^2,dims=1)))
end

"""
    rf_idx, rf_type = get_RF_types(seq, t)

Get RF centers and types (excitation or precession). Useful for k-space calculations.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) time values

# Returns
- `rf_idx`: (`::Vector{Int64}`) indices of the RF centers
- `rf_type`: (`::Vector{Int64}`, opts: [`0`, `1`]) RF type (`0`: excitation, `1`:
    precession)
"""
function get_RF_types(seq, t)
	α = get_flip_angles(seq)
	RF_mask = is_RF_on.(seq)
	RF_ex = (α .<= 90.01) .* RF_mask
	RF_rf = (α .>  90.01) .* RF_mask
	rf_idx = Int[]
	rf_type = Int[]
	T0 = cumsum([0; durs(seq)[:]])
	for i = 1:length(seq)
		if is_RF_on(seq[i])
			trf = get_RF_center(seq[i].RF[1]) + T0[i]
			append!(rf_idx, argmin(abs.(trf.-t)))
			if RF_ex[i]
				append!(rf_type, 0)
			elseif RF_rf[i]
				append!(rf_type, 1)
			end
		end
	end
	rf_idx, rf_type
end

"""
    kspace, kspace_adc = get_kspace(seq::Sequence; Δt=1)

Outputs the designed k-space trajectory of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `kspace`: (`3-column ::Matrix{Float64}`) kspace
- `kspace_adc`: (`3-column ::Matrix{Float64}`) adc kspace
"""
get_kspace(seq::Sequence; Δt=1) = begin
	t, Δt = get_uniform_times(seq, Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>[Gx; 0], 2=>[Gy; 0], 3=>[Gz; 0])
	#kspace
	Nt = length(t)
	k = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt,1,type="ordered", breaks=idx_rf)
	for i = 1:3
		kf = 0
		for (rf, p) in enumerate(parts)
			k[p,i] = cumtrapz(Δt[p]', G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1 #First part does not have RF
				k[p,i] .-= rf_type[rf-1] * kf
			end
			kf = k[p[end],i]
		end
	end
	kspace = γ * k #[m^-1]
	#Interp, as Gradients are generally piece-wise linear, the integral is piece-wise quadratic
	#Nevertheless, the integral is sampled at the ADC times so a linear interp is sufficient
	#TODO: check if this interpolation is necessary
	ts = t .+ Δt
	t_adc =  get_adc_sampling_times(seq)
	kx_adc = linear_interpolation(ts,kspace[:,1],extrapolation_bc=0)(t_adc)
	ky_adc = linear_interpolation(ts,kspace[:,2],extrapolation_bc=0)(t_adc)
	kz_adc = linear_interpolation(ts,kspace[:,3],extrapolation_bc=0)(t_adc)
	kspace_adc = [kx_adc ky_adc kz_adc]
	#Final
	kspace, kspace_adc
end

"""
    M1, M1_adc = get_M1(seq::Sequence; Δt=1)

Outputs the designed M1 of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `M1`: (`3-column ::Matrix{Float64}`) First moment
- `M1_adc`: (`3-column ::Matrix{Float64}`) First moment sampled at ADC points
"""
get_M1(seq::Sequence; Δt=1) = begin
	t, Δt = get_uniform_times(seq, Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>[Gx; 0], 2=>[Gy; 0], 3=>[Gz; 0])
	#kspace
	Nt = length(t)
	m1 = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt,1,type="ordered", breaks=idx_rf)
	for i = 1:3
		m1f = 0
		for (rf, p) in enumerate(parts)
			m1[p,i] = cumtrapz(Δt[p]', [t[p]' t[p[end]]'.+Δt[p[end]]] .* G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1 #First part does not have RF
				m1[p,i] .-= rf_type[rf-1] * m1f
			end
			m1f = m1[p[end],i]
		end
	end
	M1 = γ * m1 #[m^-1]
	#Interp, as Gradients are generally piece-wise linear, the integral is piece-wise cubic
	#Nevertheless, the integral is sampled at the ADC times so a linear interp is sufficient
	#TODO: check if this interpolation is necessary
	ts = t .+ Δt
	t_adc =  get_adc_sampling_times(seq)
	M1x_adc = linear_interpolation(ts,M1[:,1],extrapolation_bc=0)(t_adc)
	M1y_adc = linear_interpolation(ts,M1[:,2],extrapolation_bc=0)(t_adc)
	M1z_adc = linear_interpolation(ts,M1[:,3],extrapolation_bc=0)(t_adc)
	M1_adc = [M1x_adc M1y_adc M1z_adc]
	#Final
	M1, M1_adc
end


"""
    M2, M2_adc = get_M2(seq::Sequence; Δt=1)

Outputs the designed M2 of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `M2`: (`3-column ::Matrix{Float64}`) Second moment
- `M2_adc`: (`3-column ::Matrix{Float64}`) Second moment sampled at ADC points
"""
get_M2(seq::Sequence; Δt=1) = begin
	t, Δt = get_uniform_times(seq, Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>[Gx; 0], 2=>[Gy; 0], 3=>[Gz; 0])
	#kspace
	Nt = length(t)
	m2 = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt,1,type="ordered", breaks=idx_rf)
	for i = 1:3
		m2f = 0
		for (rf, p) in enumerate(parts)
			m2[p,i] = cumtrapz(Δt[p]', [t[p]' t[p[end]]'.+Δt[p[end]]].^2 .* G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1 #First part does not have RF
				m2[p,i] .-= rf_type[rf-1] * m2f
			end
			m2f = m2[p[end],i]
		end
	end
	M2 = γ * m2 #[m^-1]
	#Interp, as Gradients are generally piece-wise linear, the integral is piece-wise cubic
	#Nevertheless, the integral is sampled at the ADC times so a linear interp is sufficient
	#TODO: check if this interpolation is necessary
	ts = t .+ Δt
	t_adc =  get_adc_sampling_times(seq)
	M2x_adc = linear_interpolation(ts,M2[:,1],extrapolation_bc=0)(t_adc)
	M2y_adc = linear_interpolation(ts,M2[:,2],extrapolation_bc=0)(t_adc)
	M2z_adc = linear_interpolation(ts,M2[:,3],extrapolation_bc=0)(t_adc)
	M2_adc = [M2x_adc M2y_adc M2z_adc]
	#Final
	M2, M2_adc
end
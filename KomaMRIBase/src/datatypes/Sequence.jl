"""
    seq = Sequence()
    seq = Sequence(GR)
    seq = Sequence(GR, RF)
    seq = Sequence(GR, RF, ADC)
    seq = Sequence(GR, RF, ADC, DUR)
    seq = Sequence(GR::Array{Grad,1})
    seq = Sequence(GR::Array{Grad,1}, RF::Array{RF,1})
    seq = Sequence(GR::Array{Grad,1}, RF::Array{RF,1}, A::ADC, DUR, DEF)

The Sequence struct contains events of an MRI sequence. Most field names, except for the DEF
field, consist of matrices or vectors. In these matrices, each column index represents a
sequence block. This struct serves as input for the simulation.

# Arguments
- `GR`: (`::Matrix{Grad}`) gradient matrix. Rows for x-y-z amplitudes and columns are fo blocks
- `RF`: (`::Matrix{RF}`) RF matrix. The 1 row is for the coil and columns are for blocks
- `ADC`: (`::Array{ADC,1}`) ADC block vector
- `DUR`: (`::Vector`, `[s]`) duration block vector
- `DEF`: (`::Dict{String, Any}`) dictionary with relevant information of the sequence.
    Possible keys include [`"AdcRasterTime"`, `"GradientRasterTime"`, `"Name"`, `"Nz"`,
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
	#Ext::Array{Extension,1}
	Sequence(GR, RF, ADC, DUR, DEF) = begin
		@assert size(GR,2) .== size(RF,2) .== length(ADC) .== length(DUR) "The number of Gradient, RF, ADC, and DUR objects must be the same."
		M,N = size(GR)
		new([i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i=1:3, j=1:N],
			RF,
			ADC,
			maximum([GR.dur RF.dur ADC.dur DUR],dims=2)[:],
			DEF)
	end
end

# Main Constructors
function Sequence(GR)
	M, N = size(GR)
    gr = [i <= M ? GR[i,j] : Grad(0, 0) for i in 1:3, j in 1:N]
    rf = reshape([RF(0, 0) for i in 1:N], 1, :)
    adc = [ADC(0, 0) for _ = 1:N]
    return Sequence(gr, rf, adc, GR.dur, Dict())
end
function Sequence(GR, RF)
	@assert size(GR, 2) .== size(RF, 2) "The number of Gradient, and RF objects must be the same."
	M, N = size(GR)
    gr = [i <= M ? GR[i,j] : Grad(0, 0) for i in 1:3, j in 1:N]
    adc = [ADC(0, 0) for _ in 1:N]
    dur = maximum([GR.dur RF.dur], dims=2)[:]
	return Sequence(gr, RF, adc, dur, Dict())
end
function Sequence(GR, RF, ADC)
	@assert size(GR, 2) .== size(RF, 2) .== length(ADC) "The number of Gradient, RF, and ADC objects must be the same."
	M, N = size(GR)
    gr = [i <= M ? GR[i,j] : Grad(0, 0) for i in 1:3, j in 1:N]
    dur = maximum([GR.dur RF.dur ADC.dur], dims=2)[:]
	return Sequence(gr, RF, ADC, dur, Dict())
end
function Sequence(GR, RF, ADC, DUR)
	@assert size(GR, 2) .== size(RF, 2) .== length(ADC) .== length(DUR) "The number of Gradient, RF, ADC, and DUR objects must be the same."
	M, N = size(GR)
    gr = [i <= M ? GR[i,j] : Grad(0, GR[1,j].T, GR[1,j].rise, GR[1,j].fall, GR[1,j].delay) for i in 1:3, j in 1:N]
    dur = maximum([GR.dur RF.dur ADC.dur DUR], dims=2)[:]
	return Sequence(gr, RF, ADC, dur, Dict())
end

# Other constructors
Sequence(GR::Array{Grad,1}) = Sequence(reshape(GR,1,:))
Sequence(GR::Array{Grad,1}, RF::Array{RF,1})= Sequence(reshape(GR, :, 1), reshape(RF, 1, :), [ADC(0, 0) for i in 1:size(GR, 2)])
Sequence(GR::Array{Grad,1}, RF::Array{RF,1}, A::ADC, DUR, DEF)= Sequence(reshape(GR, :, 1), reshape(RF, 1, :), [A], [DUR], DEF)
Sequence() = Sequence([Grad(0, 0)])

# Display on the REPL
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
recursive_merge(x::AbstractDict...) = merge(recursive_merge, x...)
recursive_merge(x::AbstractVector...) = cat(x...; dims=1)
recursive_merge(x...) = x[end]
+(x::Sequence, y::Sequence) = Sequence([x.GR  y.GR], [x.RF y.RF], [x.ADC; y.ADC], [x.DUR; y.DUR], recursive_merge(x.DEF, y.DEF))
-(x::Sequence, y::Sequence) = Sequence([x.GR -y.GR], [x.RF y.RF], [x.ADC; y.ADC], [x.DUR; y.DUR], recursive_merge(x.DEF, y.DEF))
-(x::Sequence) = Sequence(-x.GR, x.RF, x.ADC, x.DUR, x.DEF)
*(x::Sequence, α::Real) = Sequence(α*x.GR, x.RF, x.ADC, x.DUR, x.DEF)
*(α::Real, x::Sequence) = Sequence(α*x.GR, x.RF, x.ADC, x.DUR, x.DEF)
*(x::Sequence, α::ComplexF64) = Sequence(x.GR, α.*x.RF, x.ADC, x.DUR, x.DEF)
*(α::ComplexF64, x::Sequence) = Sequence(x.GR, α.*x.RF, x.ADC, x.DUR, x.DEF)
*(x::Sequence, A::Matrix{Float64}) = Sequence(A*x.GR, x.RF, x.ADC, x.DUR, x.DEF) #TODO: change this, Rotation fo waveforms is broken
*(A::Matrix{Float64}, x::Sequence) = Sequence(A*x.GR, x.RF, x.ADC, x.DUR, x.DEF) #TODO: change this, Rotation fo waveforms is broken
/(x::Sequence, α::Real) = Sequence(x.GR/α, x.RF, x.ADC, x.DUR, x.DEF)
#Grad operations
+(s::Sequence, g::Grad) = s + Sequence(reshape([g],1,1)) #Changed [a;;] for reshape(a,1,1) for Julia 1.6
+(g::Grad, s::Sequence) = Sequence(reshape([g],1,1)) + s #Changed [a;;] for reshape(a,1,1) for Julia 1.6
#RF operations
+(s::Sequence, r::RF) = s + Sequence(reshape([Grad(0,0)],1,1),reshape([r],1,1)) #Changed [a;;] for reshape(a,1,1) for Julia 1.6
+(r::RF, s::Sequence) = Sequence(reshape([Grad(0,0)],1,1),reshape([r],1,1)) + s #Changed [a;;] for reshape(a,1,1) for Julia 1.6
#ADC operations
+(s::Sequence, a::ADC) = s + Sequence(reshape([Grad(0,0)],1,1),reshape([RF(0,0)],1,1),[a]) #Changed [a;;] for reshape(a,1,1) for Julia 1.6
+(a::ADC, s::Sequence) = Sequence(reshape([Grad(0,0)],1,1),reshape([RF(0,0)],1,1),[a]) + s #Changed [a;;] for reshape(a,1,1) for Julia 1.6
#Sequence object functions
size(x::Sequence) = size(x.GR[1,:])

"""
    y = is_ADC_on(x::Sequence)
    y = is_ADC_on(x::Sequence, t::AbstractVecOrMat})

Tells if the sequence `seq` has elements with ADC active, or active during time `t`.

# Arguments
- `x`: (`::Sequence`) sequence struct
- `t`: (`::AbstractVecOrMat`, `[s]`) time to check

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
    y = is_RF_on(x::Sequence, t::Vector{Real})

Tells if the sequence `seq` has elements with RF active, or active during time `t`.

# Arguments
- `x`: (`::Sequence`) Sequence struct
- `t`: (`::AbstractVector{Real}`, `[s]`) time to check

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

The total duration of the sequence in seconds.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `T`: (`::Real`, `[s]`) total duration of the sequence
"""
dur(x::Sequence) = sum(durs(x))

"""
    T0 = get_block_start_times(seq::Sequence)

Returns a vector containing the start times of blocks in a sequence. The initial time is
always zero, and the final time corresponds to the duration of the sequence.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `T0`: (`::Vector{Real}`, `[s]`) start times of the blocks in a sequence
"""
get_block_start_times(seq::Sequence) = cumsum([0; seq.DUR], dims=1)

"""
    samples = get_samples(seq::Sequence; off_val=0, max_rf_samples=Inf)

Returns the samples of the events in `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `off_val`: (`::Number`, `=0`) offset value for amplitude. Typically used to hide points in
    plots by setting it to `Inf`
- `max_rf_samples`: (`::Integer`, `=Inf`) maximum number of samples for the RF struct

# Returns
- `samples`: (`::NamedTuple`) contains samples for `gx`, `gy`, `gz`, `rf`, and `adc` events.
    Each event, represented by `e::NamedTuple`, includes time samples (`e.t`) and amplitude
    samples (`e.A`)
"""
get_samples(seq::Sequence; off_val=0, max_rf_samples=Inf) = begin
    N = length(seq)
    T0 = get_block_start_times(seq)
    # GRADs
    t_gx = reduce(vcat, [get_theo_t(seq.GR[1,i]) .+ T0[i] for i in 1:N])
    t_gy = reduce(vcat, [get_theo_t(seq.GR[2,i]) .+ T0[i] for i in 1:N])
    t_gz = reduce(vcat, [get_theo_t(seq.GR[3,i]) .+ T0[i] for i in 1:N])
    A_gx = reduce(vcat, [get_theo_A(seq.GR[1,i]; off_val) for i in 1:N])
    A_gy = reduce(vcat, [get_theo_A(seq.GR[2,i]; off_val) for i in 1:N])
    A_gz = reduce(vcat, [get_theo_A(seq.GR[3,i]; off_val) for i in 1:N])
    # RFs
    t_rf = reduce(vcat, [get_theo_t(seq.RF[1,i]; max_rf_samples) .+ T0[i] for i in 1:N])
    A_rf = reduce(vcat, [get_theo_A(rf; off_val, max_rf_samples) for rf in seq.RF])
    # ADCs
    t_adc = reduce(vcat, [get_theo_t(seq.ADC[i]) .+ T0[i] for i in 1:N])
    A_adc = reduce(vcat, [get_theo_A(adc; off_val) for adc in seq.ADC])
    return (
        gx = (t = t_gx, A = A_gx),
        gy = (t = t_gy, A = A_gy),
        gz = (t = t_gz, A = A_gz),
        rf = (t = t_rf, A = A_rf),
        adc = (t = t_adc, A = A_adc)
    )
end

"""
    y = ⏢(A, t, ΔT, ζ1, ζ2, delay)

Generates a trapezoidal waveform vector.

# Arguments
- `A`: (`::Real`) amplitude
- `t`: (`::Matrix{Real}`, `[s]`) times to evaluate (it's a `1-row matrix`)
- `ΔT`: (`::Real`, `[s]`) time duration of the top-flat
- `ζ1`: (`::Real`, `[s]`) rise time duration
- `ζ2`: (`::Real`, `[s]`) fall time duration
- `delay`: (`::Real`, `[s]`) delay time

# Returns
- `y`: (`::Matrix{Real}`) trapezoidal waveform (it's a `1-row matrix`)
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
		aux = (ζ1    .< t .- delay .< ζ1+ΔT) .* B
		if ζ1 != 0
			aux .+= (0     .< t .- delay .<= ζ1) .* A[1] .* (t .- delay)./ζ1
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

Returns the amplitude and frequency difference with respect to the resonance frequency
contained in a sequence structure.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`1-row ::Matrix{Real}`, `[s]`) time points

# Returns
- `B1`: (`1-row ::Matrix{Complex}}`, `[T]`) vector of RF pulses
- `Δf_rf`: (`1-row ::Matrix{Real}`, `[Hz]`) delta frequency vector
"""
get_rfs(seq::Sequence, t) = begin
    ϵ = MIN_RISE_TIME
    # Amplitude
    A  = seq.RF.A
    Δf = seq.RF.Δf
    # Timings
    T = seq.RF.T
    delay = seq.RF.delay
    T0 = get_block_start_times(seq)
    (sum([⏢(A[1,i], t .- T0[i], sum(T[i]) .- 2ϵ, ϵ, ϵ, delay[i]) for i=1:length(seq)]),
     sum([⏢(Δf[1,i], t .- T0[i], sum(T[i]) .- 2ϵ, ϵ, ϵ, delay[i]) for i=1:length(seq)])
    )
end

"""
    y = get_flip_angles(x::Sequence)

Returns all the flip angles of the RF pulses in the sequence `x`.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Vector{Real}`, `[deg]`) flip angles
"""
get_flip_angles(x::Sequence) = get_flip_angle.(x.RF)[:]

"""
    rf_idx, rf_type = get_RF_types(seq, t)

Get the RF centers and types (excitation or precession). Useful for k-space calculations.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Real}`, `[s]`) time values

# Returns
- `rf_idx`: (`::Vector{Integer}`) indices of the RF centers
- `rf_type`: (`::Vector{Integer}`, opts: [`0`, `1`]) RF type (`0`: excitation, `1`:
    precession)
"""
function get_RF_types(seq, t)
	α = get_flip_angles(seq)
	RF_mask = is_RF_on.(seq)
	RF_ex = (α .<= 90.01) .* RF_mask
	RF_rf = (α .>  90.01) .* RF_mask
	rf_idx = Int[]
	rf_type = Int[]
	T0 = get_block_start_times(seq)
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
- `skip_rf`: (`::Vector{Bool}`, `=zeros(Bool, sum(is_RF_on.(seq)))`) boolean vector
    indicating if the RFs should be considered for the k-space computation, with a length
    equal to the number of active RF blocks in a sequence

# Returns
- `kspace`: (`3-column ::Matrix{Real}`) k-space
- `kspace_adc`: (`3-column ::Matrix{Real}`) k-space sampled at ADC times
"""
get_kspace(seq::Sequence; Δt=1, skip_rf=zeros(Bool, sum(is_RF_on.(seq)))) = begin
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	t = t[1:end-1]
	#kspace
	Nt = length(t)
	k = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1; breaks=idx_rf)
	for i = 1:3
		kf = 0
		for (rf, p) in enumerate(parts)
			k[p,i] = cumtrapz(Δt[p]', G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1
				if !skip_rf[rf-1]
					if rf_type[rf-1] == 0 # Excitation
						k[p,i] .-= 0
					elseif rf_type[rf-1] == 1 # Refocuse
						k[p,i] .-= kf
					end
				else
					k[p,i] .+= kf
				end
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
- `skip_rf`: (`::Vector{Bool}`, `=zeros(Bool, sum(is_RF_on.(seq)))`) boolean vector
    indicating if the RFs should be considered for the M1 computation, with a length
    equal to the number of active RF blocks in a sequence

# Returns
- `M1`: (`3-column ::Matrix{Real}`) first moment
- `M1_adc`: (`3-column ::Matrix{Real}`) first moment sampled at ADC times
"""
get_M1(seq::Sequence; Δt=1, skip_rf=zeros(Bool, sum(is_RF_on.(seq)))) = begin
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	t = t[1:end-1]
	#kspace
	Nt = length(t)
	m1 = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1; breaks=idx_rf)
	for i = 1:3
		m1f = 0
		for (rf, p) in enumerate(parts)
			m1[p,i] = cumtrapz(Δt[p]', [t[p]' t[p[end]]'.+Δt[p[end]]] .* G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1 #First part does not have RF
				if !skip_rf[rf-1]
					if rf_type[rf-1] == 0 # Excitation
						m1[p,i] .-= 0
					elseif rf_type[rf-1] == 1 # Refocuse
						m1[p,i] .-= m1f
					end
				else
					m1[p,i] .+= m1f
				end
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
- `M2`: (`3-column ::Matrix{Real}`) second moment
- `M2_adc`: (`3-column ::Matrix{Real}`) second moment sampled at ADC times
"""
get_M2(seq::Sequence; Δt=1) = begin
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	t = t[1:end-1]
	#kspace
	Nt = length(t)
	m2 = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1, breaks=idx_rf)
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

"""
	SR, SR_adc = get_slew_rate(seq::Sequence; Δt=1)

Outputs the slew rate of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `SR`: (`3-column ::Matrix{Real}`) slew rate
- `SR_adc`: (`3-column ::Matrix{Real}`) slew rate sampled at ADC times
"""
get_slew_rate(seq::Sequence; Δt=1) = begin
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	t = t[1:end-1]
	#kspace
	Nt = length(t)
	m2 = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1; breaks=idx_rf)
	for i = 1:3
		m2f = 0
		for (rf, p) in enumerate(parts)
			m2[p,i] = (G[i][p[1]+1:p[end]+1] .- G[i][p[1]:p[end]])[:] ./ Δt[p]
		end
	end
	M2 = m2 #[m^-1]
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

"""
    EC, EC_adc = get_eddy_currents(seq::Sequence; Δt=1, λ=80e-3)

Outputs the Eddy currents of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients
- `λ`: (`::Real`, `=80e-3`, `[s]`) Eddy current decay constant time

# Returns
- `EC`: (`3-column ::Matrix{Real}`) Eddy currents
- `EC_adc`: (`3-column ::Matrix{Real}`) Eddy currents sampled at ADC times
"""
get_eddy_currents(seq::Sequence; Δt=1, λ=80e-3) = begin
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	t = t[1:end-1]
	#kspace
	Nt = length(t)
	m2 = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1; breaks=idx_rf)
	for i = 1:3
		m2f = 0
		for (rf, p) in enumerate(parts)
			m2[p,i] = (G[i][p[1]+1:p[end]+1] .- G[i][p[1]:p[end]])[:] ./ Δt[p]
		end
	end
	ec(t, λ) = exp.(-t ./ λ) .* (t .>= 0)
	M2 = [sum( m2[:, j] .* ec(t[i] .- t, λ) .* Δt ) for i = 1:Nt, j = 1:3] #[m^-1]
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

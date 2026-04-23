"""
    seq = Sequence()
    seq = Sequence(GR)
    seq = Sequence(GR, RF)
    seq = Sequence(GR, RF, ADC)
    seq = Sequence(GR, RF, ADC, DUR)
    seq = Sequence(GR::Array{<:Grad,1})
    seq = Sequence(GR::Array{<:Grad,1}, RF::Array{<:RF,1})
    seq = Sequence(GR::Array{<:Grad,1}, RF::Array{<:RF,1}, A::ADC, DUR, DEF)

The Sequence struct. It contains events of an MRI sequence. Most field names (except for the
DEF field) consist of matrices or vectors, where each column index represents a sequence
block. This struct serves as an input for the simulation.

# Arguments
- `GR`: (`::Matrix{Grad}`) gradient matrix. Rows for x-y-z amplitudes and columns are for blocks
- `RF`: (`::Matrix{RF}`) RF matrix. The 1 row is for the coil and columns are for blocks
- `ADC`: (`::Array{ADC,1}`) ADC block vector
- `DUR`: (`::Vector`, `[s]`) duration block vector
- `DEF`: (`::Dict{String, Any}`) dictionary with relevant information of the sequence.
    Possible keys could be [`"AdcRasterTime"`, `"GradientRasterTime"`, `"Name"`, `"Nz"`,
    `"Num_Blocks"`, `"Nx"`, `"Ny"`, `"PulseqVersion"`, `"BlockDurationRaster"`,
    `"FileName"`, `"RadiofrequencyRasterTime"`]

# Returns
- `seq`: (`::Sequence`) Sequence struct
"""
mutable struct Sequence{GT,RT,AT,DT,XT,DF}
    GR::GT #Sequence in (X, Y and Z) and time
    RF::RT #RF pulses in coil and time
    ADC::AT #ADC in time
    DUR::DT #Duration of each block, this enables delays after RF pulses to satisfy ring-down times
    EXT::XT
    DEF::DF #Dictionary with information relevant to the reconstructor
    Sequence(GR, RF, ADC, DUR, EXT, DEF) = begin
        @assert size(GR, 2) == size(RF, 2) == length(ADC) == length(DUR) "The number of Gradient, RF, ADC, DUR and EXT objects must be the same."
        GR = _ensure_gradient_axes(GR)
        new{typeof(GR),typeof(RF),typeof(ADC),typeof(DUR),typeof(EXT),typeof(DEF)}(
            GR,
            RF,
            ADC,
            DUR,
            EXT,
            DEF,
        )
    end
end

function _ensure_gradient_axes(GR)
    nchannels = size(GR, 1)
    missing_axes = 3 - nchannels
    missing_axes <= 0 && return GR
    return vcat(GR, [0.0 .* GR[1:1, :] for _ in 1:missing_axes]...)
end

_empty_extensions_per_block(nblocks) = [Extension[] for _ = 1:nblocks]

# Main Constructors
function Sequence(GR)
    rf = reshape([RF(0.0, 0.0) for i in 1:size(GR, 2)], 1, :)
    adc = [ADC(0, 0.0) for _ = 1:size(GR, 2)]
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, rf, adc, GR.dur, ext, Dict{String, Any}())
end
function Sequence(GR, RF)
    adc = [ADC(0, 0.0) for _ in 1:size(GR, 2)]
    dur = maximum([GR.dur RF.dur], dims=2)[:]
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, RF, adc, dur, ext, Dict{String, Any}())
end
function Sequence(GR, RF, ADC)
    dur = maximum([GR.dur RF.dur ADC.dur], dims=2)[:]
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, RF, ADC, dur, ext, Dict{String, Any}())
end
function Sequence(GR, RF, ADC, DUR)
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, RF, ADC, DUR, ext, Dict{String, Any}())
end
function Sequence(GR, RF, ADC, DUR, EXT)
    return Sequence(GR, RF, ADC, DUR, EXT, Dict{String, Any}())
end

# Other constructors
Sequence(GR::AbstractVector{<:Grad}) = Sequence(reshape(GR, 1, :))
Sequence(GR::AbstractVector{<:Grad}, RF::AbstractVector{<:RF}) = Sequence(reshape(GR, :, 1), reshape(RF, 1, :), [ADC(0, 0.0) for _ in 1:size(GR, 2)])
Sequence(GR::AbstractVector{<:Grad}, RF::AbstractVector{<:RF}, A::ADC, DUR, EXT, DEF) = Sequence(reshape(GR, :, 1), reshape(RF, 1, :), [A], [DUR], [EXT], DEF)
Sequence() = Sequence(
    Matrix{Grad}(undef, 3, 0),
    Matrix{RF}(undef, 1, 0),
    Vector{ADC}(undef, 0),
    Vector{Float64}(undef, 0),
    Vector{Vector{Extension}}(undef, 0),
    Dict{String, Any}(),
)

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
    if length(s) > 0
        if !compact
            nGRs = sum(is_Gx_on.(s)) + sum(is_Gy_on.(s)) + sum(is_Gz_on.(s))
            print(io, "Sequence[ τ = $(round(dur(s)*1e3;digits=3)) ms | blocks: $(length(s)) | ADC: $(sum(is_ADC_on.(s))) | GR: $nGRs | RF: $(sum(is_RF_on.(s))) | EXT: $(sum(isempty.(s.EXT) .== 0)) | DEF: $(length(s.DEF)) ]")
        else
            print(io, "Sequence[τ = $(round(dur(s)*1e3;digits=3)) ms]")
        end
    else
        print(io, "Sequence[]")
    end
end

#Sequence operations
Base.length(x::Sequence) = length(x.DUR)
Base.iterate(x::Sequence) = length(x) == 0 ? nothing : (view(x, 1), 2)
Base.iterate(x::Sequence, i::Integer) = (i <= length(x)) ? (view(x, i), i+1) : nothing
Base.view(x::Sequence, i::Integer) = (checkbounds(x.DUR, i); view(x, i:i))
Base.view(x::Sequence, i::AbstractVector{Bool}) = any(i) ? view(x, findall(i)) : nothing
Base.view(x::Sequence, i) = @views Sequence(x.GR[:,i], x.RF[:,i], x.ADC[i], x.DUR[i], x.EXT[i], x.DEF)
Base.getindex(x::Sequence, i::UnitRange{Int}) = view(x, i)
Base.getindex(x::Sequence, i::Int) = view(x, i)
Base.getindex(x::Sequence, i::AbstractVector{Bool}) = view(x, i)
Base.lastindex(x::Sequence) = length(x.DUR)
Base.copy(x::Sequence) = Sequence(_deepcopy_fields(x)...)

#Arithmetic operations
+(x::Sequence, y::Sequence) = Sequence(hcat(x.GR,  y.GR), hcat(x.RF, y.RF), vcat(x.ADC, y.ADC), vcat(x.DUR, y.DUR), vcat(x.EXT,y.EXT),merge(x.DEF, y.DEF))
-(x::Sequence, y::Sequence) = Sequence(hcat(x.GR, -y.GR), hcat(x.RF, y.RF), vcat(x.ADC, y.ADC), vcat(x.DUR, y.DUR), vcat(x.EXT,y.EXT),merge(x.DEF, y.DEF))
-(x::Sequence) = Sequence(-x.GR, x.RF, x.ADC, x.DUR, x.EXT,x.DEF)
*(x::Sequence, α::Real) = Sequence(α .* x.GR, x.RF, x.ADC, x.DUR, x.EXT, x.DEF)
*(α::Real, x::Sequence) = Sequence(α .* x.GR, x.RF, x.ADC, x.DUR, x.EXT, x.DEF)
*(x::Sequence, α::ComplexF64) = Sequence(x.GR, α.*x.RF, x.ADC, x.DUR, x.EXT, x.DEF)
*(α::ComplexF64, x::Sequence) = Sequence(x.GR, α.*x.RF, x.ADC, x.DUR, x.EXT, x.DEF)
*(x::Sequence, A::Matrix{Float64}) = Sequence(A*x.GR, x.RF, x.ADC, x.DUR, x.EXT, x.DEF) #TODO: change this, Rotation fo waveforms is broken
*(A::Matrix{Float64}, x::Sequence) = Sequence(A*x.GR, x.RF, x.ADC, x.DUR, x.EXT, x.DEF) #TODO: change this, Rotation fo waveforms is broken
/(x::Sequence, α::Real) = Sequence(x.GR/α, x.RF, x.ADC, x.DUR, x.EXT, x.DEF)
#Grad operations
+(s::Sequence, g::Grad) = s + Sequence(reshape([g],1,1)) #Changed [a;;] for reshape(a,1,1) for Julia 1.6
+(g::Grad, s::Sequence) = Sequence(reshape([g],1,1)) + s #Changed [a;;] for reshape(a,1,1) for Julia 1.6
#RF operations
+(s::Sequence, r::RF) = s + Sequence(reshape([Grad(0.0,0.0)],1,1),reshape([r],1,1)) #Changed [a;;] for reshape(a,1,1) for Julia 1.6
+(r::RF, s::Sequence) = Sequence(reshape([Grad(0.0,0.0)],1,1),reshape([r],1,1)) + s #Changed [a;;] for reshape(a,1,1) for Julia 1.6
#ADC operations
+(s::Sequence, a::ADC) = s + Sequence(reshape([Grad(0.0,0.0)],1,1),reshape([RF(0.0,0.0)],1,1),[a]) #Changed [a;;] for reshape(a,1,1) for Julia 1.6
+(a::ADC, s::Sequence) = Sequence(reshape([Grad(0.0,0.0)],1,1),reshape([RF(0.0,0.0)],1,1),[a]) + s #Changed [a;;] for reshape(a,1,1) for Julia 1.6
#Sequence object functions
size(x::Sequence) = size(x.GR[1,:])

Base.:(≈)(x::Sequence, y::Sequence; atol=1e-12) =
    typeof(x) === typeof(y) &&
    field_isapprox(x.GR, y.GR; atol=atol) &&
    field_isapprox(x.RF, y.RF; atol=atol) &&
    field_isapprox(x.ADC, y.ADC; atol=atol) &&
    field_isapprox(x.DUR, y.DUR; atol=atol) &&
    field_isapprox(x.EXT, y.EXT; atol=atol)

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
is_ADC_on(x::Sequence) = any(is_ADC_on, x.ADC)
is_ADC_on(x::Sequence, t::AbstractVecOrMat) = begin
	N = length(x)
	ts = get_block_start_times(x)[1:end-1]
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
is_RF_on(x::Sequence) = any(is_RF_on, x.RF)
is_RF_on(x::Sequence, t::AbstractVector) = begin
	N = length(x)
	ts = get_block_start_times(x)[1:end-1]
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
is_GR_on(x::Sequence) = any(is_GR_on, x.GR)

"""
    y = is_Gx_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active in x direction.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GRx in the sequence is active
"""
is_Gx_on(x::Sequence) = any(is_GR_on, x.GR.x)

"""
    y = is_Gy_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active in y direction.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GRy in the sequence is active
"""
is_Gy_on(x::Sequence) = any(is_GR_on, x.GR.y)

"""
    y = is_Gz_on(x::Sequence)

Tells if the sequence `seq` has elements with GR active in z direction.

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `y`: (`::Bool`) boolean that tells whether or not the GRz in the sequence is active
"""
is_Gz_on(x::Sequence) = any(is_GR_on, x.GR.z)

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
    T = dur(x::Sequence)

The total duration of the sequence in [s].

# Arguments
- `x`: (`::Sequence`) Sequence struct

# Returns
- `T`: (`::Real`, `[s]`) total duration of the sequence
"""
dur(x::Sequence) = sum(x.DUR)

"""
    T0 = get_block_start_times(seq::Sequence)

Returns a vector containing the start times of blocks in a sequence. The initial time is
always zero, and the final time corresponds to the duration of the sequence.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Returns
- `T0`: (`::Vector`, `[s]`) start times of the blocks in a sequence
"""
get_block_start_times(seq::Sequence) = cumsum([0.0; seq.DUR], dims=1)

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
function get_samples(seq::Sequence, range; events=[:rf, :gr, :adc], freq_in_phase=false)
    rf_samples = (;) # Empty NamedTuples
    gr_samples = (;) # Empty NamedTuples
    adc_samples = (;) # Empty NamedTuples
    T0 = get_block_start_times(seq)
    fill_if_empty(x) = isempty(x.t) && length(range) == length(seq) ? merge(x, (t=[0.0; dur(seq)], A=zeros(eltype(x.A), 2))) : x
    # RF
    if :rf in events
        t_rf = reduce(vcat, [T0[i] .+ times(seq.RF[1,i], :A)   for i in range])
        t_Δf = reduce(vcat, [T0[i] .+ times(seq.RF[1,i], :Δf)  for i in range])
        A_rf = reduce(vcat, [ampls(seq.RF[1,i]; freq_in_phase) for i in range])
        A_Δf = reduce(vcat, [freqs(seq.RF[1,i])                for i in range])
        rf_samples = (
            rf  = fill_if_empty((t = t_rf, A = A_rf)),
            Δf  = fill_if_empty((t = t_Δf, A = A_Δf))
		)
    end
    # Gradients
    if :gr in events
        t_gx = reduce(vcat, [T0[i] .+ times(seq.GR[1,i]) for i in range])
        t_gy = reduce(vcat, [T0[i] .+ times(seq.GR[2,i]) for i in range])
        t_gz = reduce(vcat, [T0[i] .+ times(seq.GR[3,i]) for i in range])
        A_gx = reduce(vcat, [ampls(seq.GR[1,i]) for i in range])
        A_gy = reduce(vcat, [ampls(seq.GR[2,i]) for i in range])
        A_gz = reduce(vcat, [ampls(seq.GR[3,i]) for i in range])
        gr_samples = (
                gx  = fill_if_empty((t = t_gx, A = A_gx)),
                gy  = fill_if_empty((t = t_gy, A = A_gy)),
                gz  = fill_if_empty((t = t_gz, A = A_gz))
                )
    end
    # ADC
    if :adc in events
        t_aq = reduce(vcat, [T0[i] .+ times(seq.ADC[i]) for i in range])
        A_aq = reduce(vcat, [ampls(seq.ADC[i]) for i in range])
        adc_samples = (
                adc = fill_if_empty((t = t_aq, A = A_aq)),
                )
    end
    # Merging events
    event_samples = merge(rf_samples, gr_samples, adc_samples)
    return event_samples
end
get_samples(seq::Sequence; kwargs...) = get_samples(seq, 1:length(seq); kwargs...)

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
function get_grads(seq, t::Union{Vector, Matrix})
    grad_samples = get_samples(seq; events=[:gr])
    for event in grad_samples
        Interpolations.deduplicate_knots!(event.t; move_knots=true)
    end
    return Tuple(linear_interpolation(event..., extrapolation_bc=0.0).(t) for event in grad_samples)
end

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
function get_rfs(seq, t::Union{Vector, Matrix})
    rf_samples = get_samples(seq; events=[:rf])
    for event in rf_samples
        Interpolations.deduplicate_knots!(event.t; move_knots=true)
    end
    return Tuple(linear_interpolation(event..., extrapolation_bc=0.0).(t) for event in rf_samples)
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
    rf_idx, rf_type = get_RF_types(seq, t)

Get RF centers and types. Useful for k-space calculations.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) time values

# Returns
- `rf_idx`: (`::Vector{Int64}`) indices of the RF centers
- `rf_types`: (`::Vector{RFUse}`) RF types
"""
function get_RF_types(seq, t)
    T0 = get_block_start_times(seq)
    rf_idx   = Int[]
    rf_types = RFUse[]
    for i in eachindex(seq.DUR)
        rf = seq.RF[1, i]
        if is_RF_on(rf)
            trf = T0[i] + rf.delay + rf.center
            push!(rf_idx, argmin(abs.(trf .- t))...)
            push!(rf_types, rf.use)
        end
    end
    return rf_idx, rf_types
end

@doc raw"""
    Mk, Mk_adc = get_Mk(seq::Sequence, k; Δt=1)

Computes the ``k``th-order moment of the Sequence `seq` given by the formula ``\int_0^T t^k G(t) dt``.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `k`: (`::Integer`) order of the moment to be computed
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `Mk`: (`3-column ::Matrix{Real}`) ``k``th-order moment
- `Mk_adc`: (`3-column ::Matrix{Real}`) ``k``th-order moment sampled at ADC times
"""
function get_Mk(seq::Sequence, k; Δt=1)
	get_sign(::Excitation) =  0
	get_sign(::Refocusing) = -1
	get_sign(::RFUse)      =  1
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	G = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	t = t[1:end-1]
	# Moment
	Nt = length(t)
	mk = zeros(Nt,3)
	# get_RF_center_breaks
	idx_rf, rf_types = get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1; breaks=idx_rf)
	for i = 1:3
		mkf = 0
		for (rf, p) in enumerate(parts)
			mk[p,i] = cumtrapz(Δt[p]', [t[p]' t[p[end]]'.+Δt[p[end]]].^k .* G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1 # First part does not have RF
				mk[p,i] .+= mkf * get_sign(rf_types[rf-1])
			end
			mkf = mk[p[end],i]
		end
	end
	Mk = γ * mk #[m^-1]
	#Interp, as Gradients are generally piece-wise linear, the integral is piece-wise cubic
	#Nevertheless, the integral is sampled at the ADC times so a linear interp is sufficient
	#TODO: check if this interpolation is necessary
	ts = t .+ Δt
	t_adc =  get_adc_sampling_times(seq)
	Mkx_adc = linear_interpolation(ts, Mk[:,1],extrapolation_bc=0)(t_adc)
	Mky_adc = linear_interpolation(ts, Mk[:,2],extrapolation_bc=0)(t_adc)
	Mkz_adc = linear_interpolation(ts, Mk[:,3],extrapolation_bc=0)(t_adc)
	Mk_adc = [Mkx_adc Mky_adc Mkz_adc]
	return Mk, Mk_adc
end

"""
Computes the k-space trajectory of the Sequence `seq`.
Refer to [`get_Mk`](@ref) and [`get_M0`](@ref)
"""
get_kspace(seq::Sequence; kwargs...) = get_Mk(seq, 0; kwargs...)

"""
Computes the zero-order moment of the Sequence `seq`.
Refer to [`get_Mk`](@ref) and [`get_kspace`](@ref)
"""
get_M0(seq::Sequence; kwargs...) = get_Mk(seq, 0; kwargs...)

"""
Computes the 1st-order moment of the Sequence `seq`.
Refer to [`get_Mk`](@ref)
"""
get_M1(seq::Sequence; kwargs...) = get_Mk(seq, 1; kwargs...)

"""
Computes the 2nd-order moment of the Sequence `seq`.
Refer to [`get_Mk`](@ref)
"""
get_M2(seq::Sequence; kwargs...) = get_Mk(seq, 2; kwargs...)

"""
	SR, SR_adc = get_slew_rate(seq::Sequence; Δt=1)

Outputs the designed slew rate of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `SR`: (`3-column ::Matrix{Float64}`) Slew rate
- `SR_adc`: (`3-column ::Matrix{Float64}`) Slew rate sampled at ADC points
"""
get_slew_rate(seq::Sequence; Δt=1) = begin
    t, Δt = get_variable_times(seq; Δt)
    Gx, Gy, Gz = get_grads(seq, t)
    t = t[1:end-1]
    Nt = length(t)
    m2 = zeros(Nt,3)
    m2[:,1] = (Gx[2:end] .- Gx[1:end-1]) ./ Δt
    m2[:,2] = (Gy[2:end] .- Gy[1:end-1]) ./ Δt
    m2[:,3] = (Gz[2:end] .- Gz[1:end-1]) ./ Δt
    Nt >= 1 && (m2[1, :] .= 0.0)
    Nt >= 2 && (m2[end, :] .= 0.0)
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

Outputs the designed eddy currents of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients
- `λ`: (`::Float64`, `=80e-3`, `[s]`) eddy current decay constant time

# Returns
- `EC`: (`3-column ::Matrix{Float64}`) Eddy currents
- `EC_adc`: (`3-column ::Matrix{Float64}`) Eddy currents sampled at ADC points
"""
get_eddy_currents(seq::Sequence; Δt=1, λ=80e-3) = begin
	t, Δt = get_variable_times(seq; Δt)
	Gx, Gy, Gz = get_grads(seq, t)
	t = t[1:end-1]
	Nt = length(t)
	m2 = zeros(Nt,3)
	m2[:,1] = (Gx[2:end] .- Gx[1:end-1]) ./ Δt
	m2[:,2] = (Gy[2:end] .- Gy[1:end-1]) ./ Δt
	m2[:,3] = (Gz[2:end] .- Gz[1:end-1]) ./ Δt
	Nt >= 1 && (m2[1, :] .= 0.0)
	Nt >= 2 && (m2[end, :] .= 0.0)
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

function get_label(seq::Sequence, nBlocks::Int=length(seq.EXT))
  if nBlocks > length(seq.EXT)
    @warn "i = $i is superior to the n° of blocks in the sequence. All the labels will be filled"
    nBlocks = length(seq.EXT)
  end

  labels = AdcLabels[]
  label = AdcLabels()

  for b = 1:nBlocks
    for val in seq.EXT[b][typeof.(seq.EXT[b]) .== LabelSet]
      setproperty!(label, Symbol(val.labelstring), val.labelvalue)
    end

    for val in seq.EXT[b][typeof.(seq.EXT[b]) .== LabelInc]
      setproperty!(label,Symbol(val.labelstring),getproperty(label,Symbol(val.labelstring)) + val.labelvalue)
    end
    push!(labels,deepcopy(label))
  end
  return labels
end

function extract_field_values(adc_labels, field)
	return [getfield(adc_label, field) for adc_label in adc_labels]
end

function Base.maximum(label::Vector{AdcLabels})
	maxLabel = AdcLabels()
	field_names = fieldnames(eltype(label))
	max_values = Dict{Symbol, Any}()

	for field in field_names
			field_values = extract_field_values(label, field)
			setproperty!(maxLabel,field,maximum(field_values))
	end

	return maxLabel
end

"""
    check_scanner_constraints(seq::Sequence, sys::Scanner=Scanner())

Checks event-local scanner limits only:
- maximum RF amplitude `|B1|`
- maximum gradient amplitude `|G|`
- maximum gradient slew rate `|dG/dt|`
- minimum ADC dwell time

It does not enforce RF/ADC dead times or RF ring-down.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `sys`: (`::Scanner`) Scanner struct
"""
function check_scanner_constraints(seq::Sequence, sys=Scanner())
    check_scanner_constraints(seq.GR, seq.RF, seq.ADC, sys)
    @info "✅ Sequence is compliant with the scanner constraints."
    return nothing
end

_absmax(x::Number) = abs(x)
_absmax(x::AbstractArray{<:Number}) = maximum(abs, x)

function check_scanner_constraints(gr::AbstractMatrix{G}, rf::AbstractMatrix{R}, adc::AbstractVector{ADC}, sys=Scanner()) where {G<:Grad,R<:RF}
    rf_maxima = IdDict{R,Float64}()
    gr_maxima = IdDict{G,Float64}()
    gr_slew_maxima = IdDict{G,Float64}()
    adc_dwells = IdDict{ADC,Float64}()
    axis_names = ("x", "y", "z")
    rtol = sqrt(eps(Float64))
    for i in eachindex(adc)
        rf_i = rf[1, i]
        if is_RF_on(rf_i) && !haskey(rf_maxima, rf_i)
            rf_maxima[rf_i] = _absmax(rf_i.A)
        end
        rf_peak = is_RF_on(rf_i) ? rf_maxima[rf_i] : 0.0
        if rf_peak > sys.B1 && !isapprox(rf_peak, sys.B1; rtol, atol=0)
            error("RF amplitude for block $i is greater than the maximum RF amplitude of the scanner ($(sys.B1 * 1e6) μT).")
        end
        for gi in 1:3
            gr_i = gr[gi, i]
            is_GR_on(gr_i) || continue
            if !haskey(gr_maxima, gr_i)
                gr_maxima[gr_i] = max(abs(gr_i.first), _absmax(gr_i.A), abs(gr_i.last))
            end
            grad_peak = gr_maxima[gr_i]
            (grad_peak <= sys.Gmax || isapprox(grad_peak, sys.Gmax; rtol, atol=0)) || error("$(axis_names[gi]) gradient amplitude for block $i is greater than the maximum gradient amplitude of the scanner ($(sys.Gmax * 1e3) mT/m).")
            if !haskey(gr_slew_maxima, gr_i)
                gr_slew_maxima[gr_i] = _max_gradient_slew(gr_i)
            end
            grad_slew = gr_slew_maxima[gr_i]
            (grad_slew <= sys.Smax || isapprox(grad_slew, sys.Smax; rtol, atol=0)) || error("$(axis_names[gi]) gradient slew rate for block $i is greater than the maximum gradient slew rate of the scanner ($(sys.Smax) mT/m/ms).")
        end
        adc_i = adc[i]
        is_ADC_on(adc_i) || continue
        if !haskey(adc_dwells, adc_i)
            adc_dwells[adc_i] = adc_i.N == 1 ? adc_i.T : adc_i.T / (adc_i.N - 1)
        end
        dwell = adc_dwells[adc_i]
        if dwell < sys.ADC_Δt
            error("ADC dwell time $(dwell * 1e6) μs for block $i is less than the minimum ADC dwell time of the scanner ($(sys.ADC_Δt * 1e6) μs).")
        end
    end
    return nothing
end

function _max_gradient_slew(gr::TrapezoidalGrad)
    max_slew = gr.rise > 0 ? abs(gr.A - gr.first) / gr.rise : 0.0
    if gr.fall > 0
        max_slew = max(max_slew, abs(gr.last - gr.A) / gr.fall)
    end
    return max_slew
end

function _max_gradient_slew(gr::UniformlySampledGrad)
    n = length(gr.A)
    n == 0 && return 0.0
    thresh = max(100 * MIN_RISE_TIME, sqrt(eps(Float64)))
    max_slew = gr.rise > 0 ? abs(gr.A[1] - gr.first) / gr.rise : 0.0
    if n > 1
        Δt = gr.T / (n - 1)
        @inbounds for i in 1:(n - 1)
            abs(Δt) > thresh || continue
            max_slew = max(max_slew, abs((gr.A[i + 1] - gr.A[i]) / Δt))
        end
    end
    if gr.fall > 0
        max_slew = max(max_slew, abs(gr.last - gr.A[end]) / gr.fall)
    end
    return max_slew
end

function _max_gradient_slew(gr::TimeShapedGrad)
    n = length(gr.A)
    n == 0 && return 0.0
    thresh = max(100 * MIN_RISE_TIME, sqrt(eps(Float64)))
    max_slew = gr.rise > 0 ? abs(gr.A[1] - gr.first) / gr.rise : 0.0
    @inbounds for i in 1:min(length(gr.T), n - 1)
        abs(gr.T[i]) > thresh || continue
        max_slew = max(max_slew, abs((gr.A[i + 1] - gr.A[i]) / gr.T[i]))
    end
    if gr.fall > 0
        max_slew = max(max_slew, abs(gr.last - gr.A[end]) / gr.fall)
    end
    return max_slew
end

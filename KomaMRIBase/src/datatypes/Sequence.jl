"""
    seq = Sequence()
    seq = Sequence(sys::Scanner)
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
- `DEF`: (`::Dict{String, Any}`) dictionary with sequence definitions. Pulseq
    raster and hardware-check metadata use Pulseq-style keys such as
    `"BlockDurationRaster"`, `"GradientRasterTime"`, `"MaxGrad"`, and `"MaxSlew"`.

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
        GR = _concrete_event_array(_ensure_gradient_axes(GR))
        RF = _concrete_event_array(RF)
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

_small_union_eltype(xs...) = _small_union_eltype_of_arrays(xs)

function _small_union_eltype_of_arrays(xs)
    types = Type[]
    eltypes = Type[]
    for x in xs
        T = eltype(x)
        push!(eltypes, T)
        if isconcretetype(T)
            T in types || push!(types, T)
        elseif T isa Union
            for U in Base.uniontypes(T)
                U in types || push!(types, U)
            end
        else
            for value in x
                T = typeof(value)
                T in types || push!(types, T)
            end
        end
    end
    return isempty(types) ? reduce(typejoin, eltypes) :
        length(types) == 1 ? only(types) : Core.apply_type(Union, types...)
end

_concrete_event_array(x) = x
_concrete_event_array(x::AbstractArray{Grad}) = _concrete_event_array_copy(x)
_concrete_event_array(x::AbstractArray{RF}) = _concrete_event_array_copy(x)

function _concrete_event_array_copy(x)
    T = _small_union_eltype(x)
    T === eltype(x) && return x
    out = similar(x, T)
    out .= x
    return out
end

const PULSEQ_RASTER_DEFINITION_KEYS = (
    "BlockDurationRaster",
    "GradientRasterTime",
    "RadiofrequencyRasterTime",
    "AdcRasterTime",
)
const PULSEQ_HW_DEFINITION_KEYS = (
    "B0",
    "MaxB1",
    "MaxGrad",
    "MaxSlew",
    "RfRingdownTime",
    "RfDeadTime",
    "AdcDeadTime",
)

const DEFAULT_SEQUENCE_DEFINITIONS = Dict{String,Any}(
    "BlockDurationRaster" => 10e-6,
    "GradientRasterTime" => 10e-6,
    "RadiofrequencyRasterTime" => 1e-6,
    "AdcRasterTime" => 100e-9,
    "B0" => 1.5,
    "MaxB1" => Inf,
    "MaxGrad" => Inf,
    "MaxSlew" => Inf,
    "RfRingdownTime" => 0.0,
    "RfDeadTime" => 0.0,
    "AdcDeadTime" => 0.0,
)

_sequence_def(sys::Scanner) = Dict{String,Any}(
    "BlockDurationRaster" => sys.DUR_Δt,
    "GradientRasterTime" => sys.GR_Δt,
    "RadiofrequencyRasterTime" => sys.RF_Δt,
    "AdcRasterTime" => sys.ADC_Δt,
    "B0" => sys.B0,
    "MaxB1" => sys.B1,
    "MaxGrad" => sys.Gmax,
    "MaxSlew" => sys.Smax,
    "RfRingdownTime" => sys.RF_ring_down_T,
    "RfDeadTime" => sys.RF_dead_time_T,
    "AdcDeadTime" => sys.ADC_dead_time_T,
)

_default_sequence_def() = copy(DEFAULT_SEQUENCE_DEFINITIONS)
_sequence_def_from_pulseq(def) = merge(_default_sequence_def(), deepcopy(def))

function _merge_sequence_def!(dest, src)
    for (key, value) in src
        haskey(dest, key) || (dest[key] = deepcopy(value))
    end
    return dest
end

function _sequence_scanner_from_def(def)
    default = DEFAULT_SEQUENCE_DEFINITIONS
    return Scanner(
        B0=get(def, "B0", default["B0"]),
        B1=get(def, "MaxB1", default["MaxB1"]),
        Gmax=get(def, "MaxGrad", default["MaxGrad"]),
        Smax=get(def, "MaxSlew", default["MaxSlew"]),
        ADC_Δt=get(def, "AdcRasterTime", default["AdcRasterTime"]),
        DUR_Δt=get(def, "BlockDurationRaster", default["BlockDurationRaster"]),
        GR_Δt=get(def, "GradientRasterTime", default["GradientRasterTime"]),
        RF_Δt=get(def, "RadiofrequencyRasterTime", default["RadiofrequencyRasterTime"]),
        RF_ring_down_T=get(def, "RfRingdownTime", default["RfRingdownTime"]),
        RF_dead_time_T=get(def, "RfDeadTime", default["RfDeadTime"]),
        ADC_dead_time_T=get(def, "AdcDeadTime", default["AdcDeadTime"]),
    )
end

# Main Constructors
function Sequence(seqs::AbstractVector{<:Sequence})
    isempty(seqs) && return Sequence()
    nblocks = sum(length, seqs)
    nblocks == 0 && return Sequence()
    i0 = findfirst(s -> length(s) > 0, seqs)
    nGR, nRF = size(seqs[i0].GR, 1), size(seqs[i0].RF, 1)
    @assert all(s -> length(s) == 0 || (size(s.GR, 1) == nGR && size(s.RF, 1) == nRF), seqs) "All sequences must have the same number of Gradient and RF channels."
    GR = Matrix{_small_union_eltype_of_arrays(s.GR for s in seqs if length(s) > 0)}(undef, nGR, nblocks)
    RF = Matrix{_small_union_eltype_of_arrays(s.RF for s in seqs if length(s) > 0)}(undef, nRF, nblocks)
    ADC = Vector{_small_union_eltype_of_arrays(s.ADC for s in seqs if length(s) > 0)}(undef, nblocks)
    DUR = Vector{Float64}(undef, nblocks)
    EXT = Vector{Vector{Extension}}(undef, nblocks)
    DEF = deepcopy(seqs[i0].DEF)
    j = 1
    for s in seqs
        n = length(s)
        n == 0 && continue
        for k in eachindex(s.DUR)
            for ax in axes(s.GR, 1)
                GR[ax, j] = copy(s.GR[ax, k])
            end
            for coil in axes(s.RF, 1)
                RF[coil, j] = copy(s.RF[coil, k])
            end
            ADC[j] = copy(s.ADC[k])
            DUR[j] = s.DUR[k]
            EXT[j] = deepcopy(s.EXT[k])
            j += 1
        end
        s === seqs[i0] || _merge_sequence_def!(DEF, s.DEF)
    end
    return Sequence(GR, RF, ADC, DUR, EXT, DEF)
end

function Sequence(GR)
    rf = reshape([RF(0.0, 0.0) for i in 1:size(GR, 2)], 1, :)
    adc = [ADC(0, 0.0) for _ = 1:size(GR, 2)]
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, rf, adc, GR.dur, ext, _default_sequence_def())
end
function Sequence(GR, RF)
    adc = [ADC(0, 0.0) for _ in 1:size(GR, 2)]
    dur = maximum([GR.dur RF.dur], dims=2)[:]
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, RF, adc, dur, ext, _default_sequence_def())
end
function Sequence(GR, RF, ADC)
    dur = maximum([GR.dur RF.dur ADC.dur], dims=2)[:]
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, RF, ADC, dur, ext, _default_sequence_def())
end
function Sequence(GR, RF, ADC, DUR)
    ext = _empty_extensions_per_block(size(GR, 2))
    return Sequence(GR, RF, ADC, DUR, ext, _default_sequence_def())
end
function Sequence(GR, RF, ADC, DUR, EXT)
    return Sequence(GR, RF, ADC, DUR, EXT, _default_sequence_def())
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
    _default_sequence_def(),
)
Sequence(sys::Scanner) = Sequence(
    Matrix{Grad}(undef, 3, 0),
    Matrix{RF}(undef, 1, 0),
    Vector{ADC}(undef, 0),
    Vector{Float64}(undef, 0),
    Vector{Vector{Extension}}(undef, 0),
    _sequence_def(sys),
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

function _copy_events(x::AbstractArray{T}) where {T}
    out = similar(x, T)
    @inbounds for i in eachindex(x)
        out[i] = copy(x[i])
    end
    return out
end
_copy_block_metadata(x) = deepcopy(x)
_fits_eltype(::Type{T}, x) where {T} = all(value -> typeof(value) <: T, x)
const _BlockEvent = Union{Grad,RF,ADC,Extension}
const _BlockEventTuple = Tuple{_BlockEvent,Vararg{_BlockEvent}}

function _event_column(events::NTuple{N,T}) where {N,T}
    out = Matrix{T}(undef, N, 1)
    @inbounds for i in 1:N
        out[i, 1] = events[i]
    end
    return out
end

function _event_column(events)
    T = _small_union_eltype(events)
    out = Matrix{T}(undef, length(events), 1)
    @inbounds for i in eachindex(events)
        out[i, 1] = events[i]
    end
    return out
end

function _event_vector(events::NTuple{N,T}) where {N,T}
    out = Vector{T}(undef, N)
    @inbounds for i in 1:N
        out[i] = events[i]
    end
    return out
end

function _event_vector(events)
    T = _small_union_eltype(events)
    out = Vector{T}(undef, length(events))
    @inbounds for i in eachindex(events)
        out[i] = events[i]
    end
    return out
end

function _append_event_matrix!(x::Matrix, y::AbstractMatrix, ::Val{copy_new}=Val(true)) where {copy_new}
    @assert size(x, 1) == size(y, 1) "Both sequences must have the same number of event channels."
    n0 = size(x, 2)
    n = size(y, 2)
    out = if _fits_eltype(eltype(x), y)
        reshape(resize!(vec(x), size(x, 1) * (n0 + n)), size(x, 1), n0 + n)
    else
        T = _small_union_eltype(x, y)
        out = Matrix{T}(undef, size(x, 1), n0 + n)
        @inbounds for j in axes(x, 2), i in axes(x, 1)
            out[i, j] = x[i, j]
        end
        out
    end
    @inbounds for j in 1:n, i in axes(y, 1)
        out[i, n0+j] = copy_new ? copy(y[i, j]) : y[i, j]
    end
    return out
end

function _append_event_vector!(x::Vector, y::AbstractVector, ::Val{copy_new}=Val(true)) where {copy_new}
    n0 = length(x)
    n = length(y)
    out = if _fits_eltype(eltype(x), y)
        resize!(x, n0 + n)
    else
        T = _small_union_eltype(x, y)
        out = Vector{T}(undef, n0 + n)
        @inbounds for j in eachindex(x)
            out[j] = x[j]
        end
        out
    end
    @inbounds for j in 1:n
        out[n0+j] = copy_new ? copy(y[j]) : y[j]
    end
    return out
end

function _append_metadata_vector!(x::Vector, y::AbstractVector, ::Val{copy_new}=Val(true)) where {copy_new}
    n0 = length(x)
    n = length(y)
    resize!(x, n0 + n)
    @inbounds for j in 1:n
        x[n0+j] = copy_new ? deepcopy(y[j]) : y[j]
    end
    return x
end

# For freshly created sequences whose events are already owned by the caller.
function _append_owned!(x::Sequence, y::Sequence)
    length(y) == 0 && return x
    x.GR = _append_event_matrix!(x.GR, y.GR, Val(false))
    x.RF = _append_event_matrix!(x.RF, y.RF, Val(false))
    x.ADC = _append_event_vector!(x.ADC, y.ADC, Val(false))
    append!(x.DUR, y.DUR)
    x.EXT = _append_metadata_vector!(x.EXT, y.EXT, Val(false))
    _merge_sequence_def!(x.DEF, y.DEF)
    return x
end

function Base.append!(x::Sequence, y::Sequence)
    length(y) == 0 && return x
    x.GR = _append_event_matrix!(x.GR, y.GR)
    x.RF = _append_event_matrix!(x.RF, y.RF)
    x.ADC = _append_event_vector!(x.ADC, y.ADC)
    append!(x.DUR, y.DUR)
    x.EXT = _append_metadata_vector!(x.EXT, y.EXT)
    _merge_sequence_def!(x.DEF, y.DEF)
    return x
end
function Base.append!(x::Sequence, y::Sequence, ys::Sequence...)
    append!(x, y)
    for z in ys
        append!(x, z)
    end
    return x
end
Base.push!(x::Sequence, y::Sequence) = append!(x, y)

_event_tuple(::Type{T}) where {T} = ()
_event_tuple(::Type{T}, event::T, rest...) where {T} = (event, _event_tuple(T, rest...)...)
_event_tuple(::Type{T}, event, rest...) where {T} = _event_tuple(T, rest...)
_check_block_event(::Union{Grad,RF,ADC,Extension}) = nothing
_check_block_event(event) = error("Unsupported block event $(typeof(event)).")
_block_duration(block_duration) = block_duration
_block_duration(block_duration, event, events...) = _block_duration(_block_duration(block_duration, event), events...)
_block_duration(block_duration, event) = max(block_duration, dur(event))
_block_duration_constraint() = nothing
_block_duration_constraint(event) = nothing
function _block_duration_constraint(event, events...)
    T = _block_duration_constraint(event)
    U = _block_duration_constraint(events...)
    if !isnothing(T) && !isnothing(U)
        error("Only one `Duration` can be added per block.")
    end
    return isnothing(T) ? U : T
end
function _set_block_duration!(seq, events...)
    block_duration = _block_duration(seq.DUR[1], events...)
    T = _block_duration_constraint(events...)
    if isnothing(T)
        seq.DUR[1] = block_duration
    elseif block_duration > T && !isapprox(block_duration, T; rtol=0, atol=PULSEQ_TIME_TOL)
        error("Block duration $(block_duration) s exceeds requested Duration($(T)) s.")
    else
        seq.DUR[1] = T
    end
    return seq
end
_axis_grad(::Nothing, ref::Nothing) = Grad(0.0, 0.0)
_axis_grad(::Nothing, ref::Grad) = 0.0 * ref
_axis_grad(gr::Grad, ref) = gr

function _block_sequence(events; x=nothing, y=nothing, z=nothing)
    foreach(_check_block_event, events)
    grads = (x, y, z)
    for (axis, gr) in zip((:x, :y, :z), grads)
        (isnothing(gr) || gr isa Grad) || error("Expected `$axis` to be a Grad.")
    end
    rfs = _event_tuple(RF, events...)
    adcs = _event_tuple(ADC, events...)
    exts = _event_tuple(Extension, events...)
    bare_grads = _event_tuple(Grad, events...)
    if !isempty(bare_grads)
        if all(isnothing, grads) && length(events) == 1 && length(bare_grads) == 1
            grads = (only(bare_grads), y, z)
        else
            error("Pass gradients with `x=`, `y=`, or `z=` when combining them with other events.")
        end
    end

    ref_gr = something(grads..., Grad(0.0, 0.0))
    grs = (_axis_grad(grads[1], ref_gr), _axis_grad(grads[2], ref_gr), _axis_grad(grads[3], ref_gr))
    rf = isempty(rfs) ? reshape([RF(0.0, 0.0)], 1, 1) : _event_column(rfs)
    adc = isempty(adcs) ? [ADC(0, 0.0)] : length(adcs) == 1 ? _event_vector(adcs) : error("Only one ADC event can be added per block.")
    seq = Sequence(_event_column(grs), rf, adc)
    _set_block_duration!(seq, events...)
    seq.EXT[1] = Extension[exts...]
    return seq
end

"""
    addblock!(seq, events...; x=nothing, y=nothing, z=nothing)

Append one block to `seq`. RF, ADC, and extensions are positional. Gradients use
`x=`, `y=`, or `z=`.
"""
function addblock!(seq::Sequence, events::_BlockEventTuple; x=nothing, y=nothing, z=nothing)
    return append!(seq, _block_sequence(events; x, y, z))
end
function addblock!(seq::Sequence, events...; x=nothing, y=nothing, z=nothing)
    return append!(seq, _block_sequence(events; x, y, z))
end
Base.push!(seq::Sequence, events...) = addblock!(seq, events...)

#Arithmetic operations
+(x::Sequence, y::Sequence) = Sequence([x, y])
-(x::Sequence, y::Sequence) = x + (-y)
-(x::Sequence) = Sequence(Matrix{Grad}(-x.GR), _copy_events(x.RF), _copy_events(x.ADC), copy(x.DUR), _copy_block_metadata(x.EXT), deepcopy(x.DEF))
*(α::Real, x::Sequence) = Sequence(Matrix{Grad}(α .* x.GR), _copy_events(x.RF), _copy_events(x.ADC), copy(x.DUR), _copy_block_metadata(x.EXT), deepcopy(x.DEF))
*(α::Complex, x::Sequence) = Sequence(_copy_events(x.GR), Matrix{RF}(α .* x.RF), Vector{ADC}(α .* x.ADC), copy(x.DUR), _copy_block_metadata(x.EXT), deepcopy(x.DEF))
function *(A::AbstractMatrix{<:Real}, x::Sequence)
    GR = Matrix{Grad}(length(x) == 1 ? reshape(A * view(x.GR, :, 1), :, 1) : A * x.GR)
    return Sequence(GR, _copy_events(x.RF), _copy_events(x.ADC), copy(x.DUR), _copy_block_metadata(x.EXT), deepcopy(x.DEF))
end
/(x::Sequence, α::Real) = Sequence(Matrix{Grad}(x.GR / α), _copy_events(x.RF), _copy_events(x.ADC), copy(x.DUR), _copy_block_metadata(x.EXT), deepcopy(x.DEF))
+(s::Sequence, events::_BlockEventTuple) = s + _block_sequence(events)
+(events::_BlockEventTuple, s::Sequence) = _block_sequence(events) + s
#Grad operations
+(s::Sequence, g::Grad) = s + (g,)
+(g::Grad, s::Sequence) = (g,) + s
#RF operations
+(s::Sequence, r::RF) = s + (r,)
+(r::RF, s::Sequence) = (r,) + s
#ADC operations
+(s::Sequence, a::ADC) = s + (a,)
+(a::ADC, s::Sequence) = (a,) + s
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

function _gradient_interpolation_samples(gr::Grad)
    t = collect(times(gr))
    A = ampls(gr)
    isempty(t) && return (t=t, A=A)
    _strictly_increasing_knots!(t)
    return (t=t, A=A)
end

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
        gx = [_gradient_interpolation_samples(seq.GR[1,i]) for i in range]
        gy = [_gradient_interpolation_samples(seq.GR[2,i]) for i in range]
        gz = [_gradient_interpolation_samples(seq.GR[3,i]) for i in range]
        t_gx = reduce(vcat, [T0[i] .+ gx[j].t for (j, i) in enumerate(range)])
        t_gy = reduce(vcat, [T0[i] .+ gy[j].t for (j, i) in enumerate(range)])
        t_gz = reduce(vcat, [T0[i] .+ gz[j].t for (j, i) in enumerate(range)])
        A_gx = reduce(vcat, [g.A for g in gx])
        A_gy = reduce(vcat, [g.A for g in gy])
        A_gz = reduce(vcat, [g.A for g in gz])
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

function get_labels(seq::Sequence, nBlocks::Int=length(seq.EXT))
  if nBlocks > length(seq.EXT)
    @warn "nBlocks = $nBlocks is greater than the number of blocks in the sequence. All labels will be filled."
    nBlocks = length(seq.EXT)
  end

  labels = AdcLabels[]
  label = AdcLabels()

  for b = 1:nBlocks
    for val in seq.EXT[b]
      label = _update_label(label, val)
    end
    push!(labels,label)
  end
  return labels
end

_update_label(label::AdcLabels, ::Extension) = label
_update_label(label::AdcLabels, val::LabelSet) = _adc_label_with(label, Symbol(val.labelstring), val.labelvalue)
function _update_label(label::AdcLabels, val::LabelInc)
  label_symbol = Symbol(val.labelstring)
  return _adc_label_with(label, label_symbol, getproperty(label, label_symbol) + val.labelvalue)
end

function Base.maximum(label::Vector{AdcLabels})
	isempty(label) && throw(ArgumentError("collection must be non-empty"))
	values = ntuple(i -> maximum(getfield(adc_label, i) for adc_label in label), Val(length(ADC_LABEL_NAMES)))
	return AdcLabels(values...)
end

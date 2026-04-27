
"""
    collect_pulseq_assets(seq::Sequence, raster::PulseqRaster) -> (blocks, event_libraries)

Create the Pulseq export dictionaries required to serialize `seq` into the Pulseq file format.
This function is responsible for deduplicating reusable objects (RF, gradients, shapes, etc.)
and for translating each sequence block into integer lookups expected by the specification.
"""
struct PulseqShapeCacheEntry
    shape_id::Int
    quantized::Vector{Int}
end

struct PulseqTrapGradEventKey
    amplitude::Float64
    rise::Float64
    flat::Float64
    fall::Float64
    delay::Float64
end

struct PulseqArbGradEventKey
    amplitude::Float64
    first::Float64
    last::Float64
    amp_shape_id::Int
    time_shape_id::Int
    delay::Float64
end

# Dedup tables for emitted event ids and reusable shape payloads.
struct PulseqWriteCache{G<:Grad,R<:RF}
    rf_event_ids::Dict{PulseqRFEvent,Int}
    trap_grad_event_ids::Dict{PulseqTrapGradEventKey,Int}
    arb_grad_event_ids::Dict{PulseqArbGradEventKey,Int}
    adc_event_ids::Dict{PulseqADCEvent,Int}
    shape_ids::Dict{Tuple{Int,UInt},Vector{PulseqShapeCacheEntry}}
    phase_shape_ids::Dict{Tuple{Int,UInt},Vector{PulseqShapeCacheEntry}}
    rf_magnitude_shape_ids::Dict{Int,Vector{PulseqShapeCacheEntry}}
    rf_object_ids::Dict{R,Int}
    grad_object_ids::Dict{G,Int}
    adc_object_ids::Dict{ADC,Int}
end

PulseqWriteCache(::Type{G}, ::Type{R}) where {G<:Grad,R<:RF} = PulseqWriteCache(
    Dict{PulseqRFEvent,Int}(),
    Dict{PulseqTrapGradEventKey,Int}(),
    Dict{PulseqArbGradEventKey,Int}(),
    Dict{PulseqADCEvent,Int}(),
    Dict{Tuple{Int,UInt},Vector{PulseqShapeCacheEntry}}(),
    Dict{Tuple{Int,UInt},Vector{PulseqShapeCacheEntry}}(),
    Dict{Int,Vector{PulseqShapeCacheEntry}}(),
    Dict{R,Int}(),
    Dict{G,Int}(),
    Dict{ADC,Int}(),
)

const PULSEQ_NO_EVENT_ID = 0
const PULSEQ_DEFAULT_TIME_SHAPE_ID = 0
const PULSEQ_OVERSAMPLED_TIME_SHAPE_ID = -1
# Dedup keys ignore roundoff; stored event values stay unchanged.
const PULSEQ_EVENT_SIGNIFICANT_DIGITS = 12

_pulseq_absmax(x::Number) = abs(x)
_pulseq_absmax(x::AbstractArray{<:Number}) = maximum(abs, x)
_is_pulseq_rf_on(rf) = _pulseq_absmax(rf.A) > PULSEQ_SHAPE_ZERO_TOL
function _is_pulseq_grad_on(grad)
    peak = max(abs(grad.first), _pulseq_absmax(grad.A), abs(grad.last))
    return iszero(peak) ? !iszero(dur(grad)) : peak > PULSEQ_SHAPE_ZERO_TOL
end

collect_pulseq_assets(seq, raster) =
    collect_pulseq_assets(seq.GR, seq.RF, seq.ADC, seq.DUR, seq.EXT, seq.DEF, raster)

pulseq_raster_tuple(raster) = (;
    BlockDurationRaster=raster.BlockDurationRaster,
    GradientRasterTime=raster.GradientRasterTime,
    RadiofrequencyRasterTime=raster.RadiofrequencyRasterTime,
    AdcRasterTime=raster.AdcRasterTime,
)

function pulseq_data(seq::KomaMRIBase.Sequence, raster::PulseqRaster)
    prepared = prepare_pulseq_write(seq, raster)
    blocks, event_libraries = collect_pulseq_assets(prepared, raster)
    return PulseqSequenceData(blocks, event_libraries, v"1.5.1", nothing)
end

function collect_pulseq_assets(gr, rf, adc, block_durations, ext, def, raster)
    blocks = Vector{PulseqBlockEventIDs}(undef, length(block_durations))
    rf_library = Dict{Int,PulseqRFEvent}()
    grad_library = Dict{Int,PulseqGradEvent}()
    adc_library = Dict{Int,PulseqADCEvent}()
    tmp_delay_library = Dict{Int,Float64}()
    shape_library = ShapeLibrary()
    extension_instance_library = Dict{Int,PulseqExtensionInstanceEvent}()
    extension_type_library = Dict{Int,Type{<:Extension}}()
    extension_spec_library = Dict{Int,Dict{Int,Extension}}()
    cache = PulseqWriteCache(eltype(gr), eltype(rf))
    extension_vector_ids = Dict{Tuple{Vararg{Extension}},Int}()
    for idx in eachindex(block_durations)
        rf_i = rf[1, idx]
        rf_id = register_rf!(rf_library, shape_library, cache, rf_i, raster)

        gx = gr[1, idx]
        gy = gr[2, idx]
        gz = gr[3, idx]
        gx_id = register_grad!(grad_library, shape_library, cache, gx, raster)
        gy_id = register_grad!(grad_library, shape_library, cache, gy, raster)
        gz_id = register_grad!(grad_library, shape_library, cache, gz, raster)

        adc_i = adc[idx]
        adc_id = register_adc!(adc_library, cache, adc_i)
        duration = pulseq_block_duration_ticks(block_durations[idx], raster.BlockDurationRaster)

        ext_vec = ext[idx]
        ext_id = isempty(ext_vec) ? PULSEQ_NO_EVENT_ID : get!(extension_vector_ids, Tuple(ext_vec)) do
            register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext_vec)
        end
        blocks[idx] = PulseqBlockEventIDs(idx, duration, rf_id, gx_id, gy_id, gz_id, adc_id, ext_id)
    end

    definitions = pulseq_definitions(def, raster)
    for ext_type in values(extension_type_library)
        ext_type === QuaternionRot || continue
        "ROTATIONS" in definitions.required_extensions || push!(definitions.required_extensions, "ROTATIONS")
    end
    event_libraries = PulseqFileEventLibraries(
        grad_library,
        rf_library,
        adc_library,
        tmp_delay_library,
        shape_library,
        extension_instance_library,
        extension_type_library,
        extension_spec_library,
        definitions,
    )
    return blocks, event_libraries
end

"""
    id = register_rf!(rf_library, shape_library, rf, raster)
"""
function register_rf!(rf_library, shape_library, cache, rf, raster)
    _is_pulseq_rf_on(rf) || return PULSEQ_NO_EVENT_ID
    event_id = get(cache.rf_object_ids, rf, PULSEQ_NO_EVENT_ID)
    event_id != PULSEQ_NO_EVENT_ID && return event_id
    rf_raster_time = raster.RadiofrequencyRasterTime
    mag_id, phase_id, time_id = register_rf_shapes!(shape_library, cache, rf, rf_raster_time)
    first_sample_offset = pulseq_rf_first_sample_offset(time_id, rf_raster_time)
    delay = rf.delay - first_sample_offset
    amp = pulseq_rf_amplitude(rf)
    freq = pulseq_rf_frequency(rf)
    phase = rf.ϕ
    use_char = get_char_from_RF_use(rf.use)
    center = isnothing(rf.center) ? nothing : rf.center + first_sample_offset
    rf_event = PulseqRFEvent(amp, mag_id, phase_id, time_id, center, delay, 0.0, 0.0, freq, phase, use_char)
    event_id = _store_event_cached!(rf_library, cache.rf_event_ids, rf_event)
    cache.rf_object_ids[rf] = event_id
    return event_id
end

pulseq_rf_amplitude(rf::BlockPulseRF) = abs(rf.A)
pulseq_rf_amplitude(rf::RF) = maximum(abs, rf.A)
pulseq_rf_frequency(rf::RF) = rf.Δf
pulseq_rf_frequency(::FrequencyModulatedRF) =
    throw(ArgumentError("Pulseq write does not support RF frequency waveforms."))

pulseq_rf_phase_shape(::Number) = [0.0, 0.0]
pulseq_rf_phase_shape(A::AbstractVector) = map(A) do a
    iszero(a) ? 0.0 : mod(angle(a), 2π) / 2π
end

pulseq_rf_magnitude_shape(::Number) = [1.0, 1.0]
function pulseq_rf_magnitude_shape(A::AbstractVector)
    mags = abs.(A)
    peak = maximum(mags)
    return iszero(peak) ? zero.(mags) : mags ./ peak
end

function register_rf_magnitude_shape!(shapes, cache::PulseqWriteCache, A)
    samples = pulseq_rf_magnitude_shape(A)
    quantized = quantized_shape_samples(samples)
    candidates = get(cache.rf_magnitude_shape_ids, length(samples), nothing)
    if !isnothing(candidates)
        shape_id = matching_rf_magnitude_shape_id(candidates, quantized)
        !isnothing(shape_id) && return shape_id
    end
    shape_id = _store_shape!(shapes, cache, samples)
    entries = get!(cache.rf_magnitude_shape_ids, length(samples), PulseqShapeCacheEntry[])
    any(entry -> entry.shape_id == shape_id, entries) || push!(entries, PulseqShapeCacheEntry(shape_id, quantized))
    return shape_id
end

# MATLAB Pulseq detects the default RF timing by checking whether the compressed
# `event.t / rfRasterTime` matches the compact center-sampled pattern in
# `matlab/+mr/@Sequence/Sequence.m` (v1.5.1, lines 562-567). In Koma we store RF
# delay to the first sample, so the same compact encoding is recognized from the
# first-sample-relative times together with the half-raster delay shift below.
function pulseq_rf_compact_time_shape_id(::TimeShapedRF, rf_raster_time)
    return nothing
end

function pulseq_rf_compact_time_shape_id(rf, rf_raster_time)
    rf.delay + PULSEQ_TIME_TOL < rf_raster_time / 2 && return nothing
    tt = _shape_times(rf.A, rf.T)
    default_tt = (0:(length(tt) - 1)) .* rf_raster_time
    if all(isapprox.(tt, default_tt; rtol=0, atol=PULSEQ_TIME_TOL))
        return PULSEQ_DEFAULT_TIME_SHAPE_ID
    end
    oversampled_tt = (0:(length(tt) - 1)) .* (rf_raster_time / 2)
    if isodd(length(tt)) && all(isapprox.(tt, oversampled_tt; rtol=0, atol=PULSEQ_TIME_TOL))
        return PULSEQ_OVERSAMPLED_TIME_SHAPE_ID
    end
    return nothing
end

function pulseq_rf_first_sample_offset(time_shape_id, rf_raster_time)
    return time_shape_id <= PULSEQ_DEFAULT_TIME_SHAPE_ID ? rf_raster_time / 2 : 0.0
end

function pulseq_rf_timing_policy(rf, Δt_rf)
    time_shape_id = pulseq_rf_compact_time_shape_id(rf, Δt_rf)
    time_step = time_shape_id == PULSEQ_OVERSAMPLED_TIME_SHAPE_ID ? Δt_rf / 2 : Δt_rf
    first_sample_offset = isnothing(time_shape_id) ? 0.0 : pulseq_rf_first_sample_offset(time_shape_id, Δt_rf)
    return (; time_shape_id, time_step, first_sample_offset)
end

function pulseq_rf_time_shape_id!(shape_library, cache, rf, rf_raster_time)
    policy = pulseq_rf_timing_policy(rf, rf_raster_time)
    if !isnothing(policy.time_shape_id)
        pulseq_delay = rf.delay - policy.first_sample_offset
        rounded_delay = round(pulseq_delay / rf_raster_time) * rf_raster_time
        isapprox(pulseq_delay, rounded_delay; rtol=0, atol=PULSEQ_TIME_TOL) && return policy.time_shape_id
    end
    return _store_shape!(shape_library, cache, _shape_times(rf.A, rf.T) ./ rf_raster_time)
end

function register_rf_shapes!(shapes, cache, rf, rf_raster_time)
    mag_id      = register_rf_magnitude_shape!(shapes, cache, rf.A)
    phase_id    = _store_shape!(shapes, cache, pulseq_rf_phase_shape(rf.A); phase_shape=true)
    time_id     = pulseq_rf_time_shape_id!(shapes, cache, rf, rf_raster_time)
    return mag_id, phase_id, time_id
end

"""
    id = register_grad!(grad_library, shape_library, grad, raster)
"""
function register_grad!(grad_library, shape_library, cache, grad, raster)
    _is_pulseq_grad_on(grad) || return PULSEQ_NO_EVENT_ID
    event_id = get(cache.grad_object_ids, grad, PULSEQ_NO_EVENT_ID)
    event_id != PULSEQ_NO_EVENT_ID && return event_id
    event_id = register_grad_event!(grad_library, shape_library, cache, grad, raster)
    cache.grad_object_ids[grad] = event_id
    return event_id
end

function register_grad_event!(grad_library, shape_library, cache, grad::TrapezoidalGrad, raster)
    if iszero(grad.first) && iszero(grad.last)
        return _store_grad_event!(grad_library, cache, PulseqTrapGradEvent(grad.A, grad.rise, grad.T, grad.fall, grad.delay))
    end
    return register_grad_event!(grad_library, shape_library, cache, edge_timed_grad(grad), raster)
end

function register_grad_event!(grad_library, shape_library, cache, grad::UniformlySampledGrad, raster)
    trap_event = pulseq_trap_event(grad)
    !isnothing(trap_event) && return _store_grad_event!(grad_library, cache, trap_event)
    Δt_gr = raster.GradientRasterTime
    amp = pulseq_gradient_amplitude(grad.A)
    shape = iszero(amp) ? zero.(grad.A) : grad.A ./ amp
    shape_id = _store_shape!(shape_library, cache, shape)
    time_id = pulseq_gradient_time_shape_id!(shape_library, cache, pulseq_grad_sample_times(grad), Δt_gr)
    return _store_grad_event!(grad_library, cache, PulseqArbGradEvent(amp, grad.first, grad.last, shape_id, time_id, grad.delay))
end

function register_grad_event!(grad_library, shape_library, cache, grad::TimeShapedGrad, raster)
    trap_event = pulseq_trap_event(grad)
    !isnothing(trap_event) && return _store_grad_event!(grad_library, cache, trap_event)
    amp = pulseq_gradient_amplitude(grad.A)
    shape = iszero(amp) ? zero.(grad.A) : grad.A ./ amp
    shape_id = _store_shape!(shape_library, cache, shape)
    time_id = pulseq_gradient_time_shape_id!(shape_library, cache, pulseq_grad_sample_times(grad), raster.GradientRasterTime)
    return _store_grad_event!(grad_library, cache, PulseqArbGradEvent(amp, grad.first, grad.last, shape_id, time_id, grad.delay))
end

edge_timed_grad(grad::TrapezoidalGrad) =
    Grad([grad.first, grad.A, grad.A, grad.last], [grad.rise, grad.T, grad.fall], 0.0, 0.0, grad.delay, grad.first, grad.last)

function pulseq_grad_sample_times(grad::UniformlySampledGrad)
    n = length(grad.A)
    interval = n > 1 ? grad.T / (n - 1) : 0.0
    return grad.rise .+ (0:(n - 1)) .* interval
end

pulseq_grad_sample_times(grad::TimeShapedGrad) =
    grad.rise .+ cumsum(vcat(0.0, grad.T))

function pulseq_trap_event(grad::UniformlySampledGrad)
    iszero(grad.first) && iszero(grad.last) || return nothing
    length(grad.A) == 2 || return nothing
    grad.A[1] == grad.A[2] || return nothing
    return PulseqTrapGradEvent(grad.A[1], grad.rise, grad.T, grad.fall, grad.delay)
end

function pulseq_trap_event(grad::TimeShapedGrad)
    iszero(grad.first) && iszero(grad.last) || return nothing
    length(grad.A) == 2 || return nothing
    length(grad.T) == 1 || return nothing
    grad.A[1] == grad.A[2] || return nothing
    return PulseqTrapGradEvent(grad.A[1], grad.rise, grad.T[1], grad.fall, grad.delay)
end

pulseq_trap_event(::Grad) = nothing

"""
    id = register_adc!(adc_library, adc)

In Pulseq the ADC event starts at a dwell-time edge, but the first sample is taken at dwell/2
(see [Pulseq time and shape specification](https://pulseq.github.io/pulseq_shapes_and_times.pdf)).
Koma uses delay = time to first sample, so we store pulseq_delay = delay - dwell_s/2.
"""
function register_adc!(adc_library, cache, adc)
    is_on(adc) || return PULSEQ_NO_EVENT_ID
    event_id = get(cache.adc_object_ids, adc, PULSEQ_NO_EVENT_ID)
    event_id != PULSEQ_NO_EVENT_ID && return event_id
    dwell_s = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
    delay_s = adc.delay - dwell_s / 2
    event = PulseqADCEvent(adc.N, dwell_s, delay_s, 0.0, 0.0, adc.Δf, adc.ϕ, 0)
    event_id = _store_event_cached!(adc_library, cache.adc_event_ids, event)
    cache.adc_object_ids[adc] = event_id
    return event_id
end

"""
    id = register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext)
"""
function register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext)
    isempty(ext) && return PULSEQ_NO_EVENT_ID
    instance_ids = Vector{Int}(undef, length(ext))
    for (i, e) in pairs(ext)
        ext_id = _store_event!(extension_type_library, typeof(e))
        if !haskey(extension_spec_library, ext_id)
            extension_spec_library[ext_id] = Dict{Int,Extension}()
        end
        ref = _store_event!(extension_spec_library[ext_id], e)
        instance_ids[i] = _store_event!(extension_instance_library, PulseqExtensionInstanceEvent((ext_id, ref, PULSEQ_NO_EVENT_ID)))
    end

    for i in firstindex(instance_ids):(lastindex(instance_ids) - 1)
        current_id = instance_ids[i]
        next_id = instance_ids[i + 1]
        current = extension_instance_library[current_id]
        extension_instance_library[current_id] = PulseqExtensionInstanceEvent(current.type, current.ref, next_id)
    end

    return instance_ids[1]
end

# store_event! and store_shape! are helper functions to store events and shapes in the dictionary, and return the id of the event or shape.
function _store_event!(event_dict::AbstractDict, event)
    for (k, v) in event_dict
        if v == event return k end
    end
    new_id = length(event_dict) + 1
    event_dict[new_id] = event
    return new_id
end

_pulseq_event_float(x) = round(x; sigdigits=PULSEQ_EVENT_SIGNIFICANT_DIGITS)
_pulseq_event_edge(x, amp) = abs(x) <= abs(amp) * PULSEQ_SHAPE_ZERO_TOL ? 0.0 : _pulseq_event_float(x)

_pulseq_grad_event_key(event::PulseqTrapGradEvent) = PulseqTrapGradEventKey(
    _pulseq_event_float(event.amplitude),
    event.rise,
    event.flat,
    event.fall,
    event.delay,
)

_pulseq_grad_event_key(event::PulseqArbGradEvent) = PulseqArbGradEventKey(
    _pulseq_event_float(event.amplitude),
    _pulseq_event_edge(event.first, event.amplitude),
    _pulseq_event_edge(event.last, event.amplitude),
    event.amp_shape_id,
    event.time_shape_id,
    event.delay,
)

function _store_event_cached!(event_dict::AbstractDict, event_cache::AbstractDict, key, event=key)
    get!(event_cache, key) do
        new_id = length(event_dict) + 1
        event_dict[new_id] = event
        new_id
    end
end

_canonical_pulseq_grad_event(event::PulseqTrapGradEvent) = event
_canonical_pulseq_grad_event(event::PulseqArbGradEvent) = PulseqArbGradEvent(
    event.amplitude,
    _pulseq_event_edge(event.first, event.amplitude),
    _pulseq_event_edge(event.last, event.amplitude),
    event.amp_shape_id,
    event.time_shape_id,
    event.delay,
)

_store_grad_event!(grad_library, cache::PulseqWriteCache, event::PulseqTrapGradEvent) =
    _store_event_cached!(grad_library, cache.trap_grad_event_ids, _pulseq_grad_event_key(event), event)
_store_grad_event!(grad_library, cache::PulseqWriteCache, event::PulseqArbGradEvent) =
    _store_event_cached!(grad_library, cache.arb_grad_event_ids, _pulseq_grad_event_key(event), _canonical_pulseq_grad_event(event))

function _store_shape!(
    shapes,
    shape_cache,
    samples;
    phase_shape = false
)
    key = shape_cache_key(samples; phase_shape)
    candidates = get(shape_cache, key, nothing)
    if !isnothing(candidates)
        shape_id = matching_shape_id(candidates, samples; phase_shape)
        !isnothing(shape_id) && return shape_id
    end

    payload = compress_shape(samples)
    new_id = length(shapes) + 1
    shapes[new_id] = payload
    push!(get!(shape_cache, key, PulseqShapeCacheEntry[]), PulseqShapeCacheEntry(new_id, quantized_shape_samples(samples; phase_shape)))
    return new_id
end

_store_shape!(shapes, cache::PulseqWriteCache, samples; phase_shape = false) =
    _store_shape!(
        shapes,
        phase_shape ? cache.phase_shape_ids : cache.shape_ids,
        samples;
        phase_shape,
    )

function get_quantized!(builder, cache, object)
    haskey(cache, object) && return cache[object]
    value = builder()
    cache[object] = value
    return value
end

function canonicalize_quantized!(objects::AbstractArray{T}, cache::IdDict{T,T}) where {T}
    for i in eachindex(objects)
        object = objects[i]
        objects[i] = get!(cache, object, object)
    end
    return objects
end

shape_cache_key(samples; phase_shape = false) =
    (length(samples), shape_cache_hash(samples; phase_shape))

function shape_cache_hash(samples; phase_shape = false)
    offset = phase_shape ? first(samples) : 0.0
    h = UInt(0)
    for sample in samples
        h = hash(quantized_shape_sample(sample, offset; phase_shape), h)
    end
    return h
end

function matching_shape_id(candidates, samples; phase_shape = false)
    offset = phase_shape ? first(samples) : 0.0
    for candidate in candidates
        length(candidate.quantized) == length(samples) || continue
        matches = true
        for i in eachindex(candidate.quantized, samples)
            if candidate.quantized[i] != quantized_shape_sample(samples[i], offset; phase_shape)
                matches = false
                break
            end
        end
        matches && return candidate.shape_id
    end
    return nothing
end

function quantized_shape_samples(samples; phase_shape = false)
    offset = phase_shape ? first(samples) : 0.0
    quantized = Vector{Int}(undef, length(samples))
    for i in eachindex(samples)
        quantized[i] = quantized_shape_sample(samples[i], offset; phase_shape)
    end
    return quantized
end

quantized_shape_sample(sample, offset; phase_shape = false) =
    round(Int, (phase_shape ? mod(sample - offset, 1.0) : sample) / PULSEQ_SHAPE_QUANTIZATION)

function matching_rf_magnitude_shape_id(candidates, quantized)
    for candidate in candidates
        length(candidate.quantized) == length(quantized) || continue
        all(abs(candidate.quantized[i] - quantized[i]) <= 1 for i in eachindex(quantized)) || continue
        return candidate.shape_id
    end
    return nothing
end

function pulseq_block_duration_ticks(duration, raster)
    ticks = duration / raster
    rounded = round(Int, ticks)
    rounded_duration = rounded * raster
    return isapprox(duration, rounded_duration; rtol=0, atol=PULSEQ_TIME_TOL) ? rounded : ceil(Int, ticks - PULSEQ_TIME_TOL / raster)
end

function pulseq_gradient_amplitude(A)
    amp = maximum(abs, A)
    iszero(amp) && return amp
    idx = findfirst(a -> abs(a) > amp * PULSEQ_SHAPE_QUANTIZATION / 2, A)
    return amp * sign(A[something(idx, findfirst(!iszero, A))])
end

function pulseq_gradient_time_shape_id!(shape_library, cache, tt, Δt)
    n = length(tt)
    n == 0 && return PULSEQ_DEFAULT_TIME_SHAPE_ID
    default_tt = ((1:n) .- 0.5) .* Δt
    if all(isapprox.(tt, default_tt; rtol=0, atol=PULSEQ_TIME_TOL))
        return PULSEQ_DEFAULT_TIME_SHAPE_ID
    end
    oversampled_tt = (1:n) .* (Δt / 2)
    if isodd(n) && all(isapprox.(tt, oversampled_tt; rtol=0, atol=PULSEQ_TIME_TOL))
        return PULSEQ_OVERSAMPLED_TIME_SHAPE_ID
    end
    return _store_shape!(shape_library, cache, tt ./ Δt)
end

"""
    emit_pulseq(io::IO, blocks, event_libraries)

Write the Pulseq sections to `io`, using the already prepared `blocks` and `event_libraries`.
"""
function emit_pulseq(io, blocks, event_libraries)
    has_arb_grads = false
    has_trap_grads = false
    for grad in values(event_libraries.grad_library)
        has_arb_grads |= grad isa PulseqArbGradEvent
        has_trap_grads |= grad isa PulseqTrapGradEvent
        has_arb_grads && has_trap_grads && break
    end
    emit_header_comment!(io)
    emit_version_section!(io)
    emit_definitions_section!(io, event_libraries.definitions)
    if isempty(blocks)
        @warn "Pulseq export: block table is empty; populate it via `collect_pulseq_assets`."
    else
        emit_blocks_section!(io, blocks)
    end
    if !isempty(event_libraries.rf_library)
        emit_rf_section!(io, event_libraries)
    end
    if has_arb_grads
        emit_gradients_section!(io, event_libraries)
    end
    if has_trap_grads
        emit_trap_section!(io, event_libraries)
    end
    if !isempty(event_libraries.adc_library)
        emit_adc_section!(io, event_libraries)
    end
    if !isempty(event_libraries.extension_instance_library)
        emit_extension_section!(io, event_libraries)
    end
    if !isempty(event_libraries.shape_library)
        emit_shapes_section!(io, event_libraries)
    end
end

emit_pulseq(io, data::PulseqSequenceData) =
    emit_pulseq(io, data.blocks, data.libraries)

function emit_header_comment!(io)
    write(io, "# Pulseq sequence file\n")
    write(io, "# Created by KomaMRI.jl\n")
    write(io, "# KomaMRIFiles v$(pkgversion(KomaMRIFiles))\n\n")
    return nothing
end

function emit_version_section!(io, version=v"1.5.1")
    write(io, "[VERSION]\n")
    write(io, "major $(version.major)\n")
    write(io, "minor $(version.minor)\n")
    write(io, "revision $(version.patch)\n\n")
end

dense_event_ids(event_dict) = Base.OneTo(length(event_dict))
emit_pulseq_float!(io, value::Real) = @printf(io, "%.12g", value)
pulseq_float_string(value::Real) = @sprintf("%.12g", value)

function emit_definitions_section!(io::IO, definitions)
    write(io, "[DEFINITIONS]\n")
    for (key, value) in definitions.entries
        should_emit_definition(key, value) || continue
        write(io, key)
        write(io, " ")
        write(io, value)
        write(io, "\n")
    end
    write(io, "BlockDurationRaster "); emit_pulseq_float!(io, definitions.block_duration_raster); write(io, "\n")
    write(io, "GradientRasterTime "); emit_pulseq_float!(io, definitions.gradient_raster_time); write(io, "\n")
    write(io, "RadiofrequencyRasterTime "); emit_pulseq_float!(io, definitions.radiofrequency_raster_time); write(io, "\n")
    write(io, "AdcRasterTime "); emit_pulseq_float!(io, definitions.adc_raster_time); write(io, "\n")
    isempty(definitions.required_extensions) || begin
        write(io, "RequiredExtensions ")
        emit_definition_value!(io, definitions.required_extensions)
        write(io, "\n")
    end
    write(io, "\n")
end

function emit_blocks_section!(io, blocks)
    write(io, "# Format of blocks:\n")
    write(io, "# NUM DUR RF  GX  GY  GZ  ADC  EXT\n")
    write(io, "[BLOCKS]\n")
    isempty(blocks) && return
    id_width = ndigits(length(blocks))
    for block in blocks
        @printf(io, "%*d %3d %3d %3d %3d %3d %2d %2d\n", id_width, block.block_id, block.duration_ticks, block.rf_id, block.gx_id, block.gy_id, block.gz_id, block.adc_id, block.ext_id)
    end
end

function emit_rf_section!(io::IO, event_libraries)
    write(io, "\n# Format of RF events:\n")
    write(io, "# id amp mag_id phase_id time_id center delay freq_ppm phase_ppm freq_off phase_off use\n")
    write(io, "# ..  Hz     ..       ..       ..    us    us      ppm   rad/MHz       Hz       rad  ..\n")
    write(io, "# Field 'use' is the initial of:\n")
    write(io, "# excitation refocusing inversion saturation preparation other undefined\n")
    write(io, "[RF]\n")
    isempty(event_libraries.rf_library) && return
    for id in dense_event_ids(event_libraries.rf_library)
        rf_data = event_libraries.rf_library[id]
        center_us = pulseq_rf_center_for_write(rf_data, event_libraries) * 1e6
        delay_us = round(Int, rf_data.delay * 1e6)
        @printf(io, "%d ", id)
        emit_pulseq_float!(io, γ * rf_data.amplitude)
        @printf(io, " %d %d %d ", rf_data.mag_id, rf_data.phase_id, rf_data.time_shape_id)
        emit_pulseq_float!(io, center_us)
        @printf(io, " %d ", delay_us)
        emit_pulseq_float!(io, rf_data.freq_ppm); write(io, " ")
        emit_pulseq_float!(io, rf_data.phase_ppm); write(io, " ")
        emit_pulseq_float!(io, rf_data.freq); write(io, " ")
        emit_pulseq_float!(io, rf_data.phase)
        @printf(io, " %c\n", rf_data.use)
    end
end

function pulseq_rf_center_for_write(rf_data, event_libraries)
    !isnothing(rf_data.center) && return rf_data.center
    Δt_rf = event_libraries.definitions.radiofrequency_raster_time
    rf = get_RF(rf_data, event_libraries.shape_library, Δt_rf)
    return rf.center + pulseq_rf_first_sample_offset(rf_data, event_libraries.shape_library, Δt_rf)
end

function pulseq_rf_first_sample_offset(rf_data, shape_library, Δt_rf)
    time_shape_id = rf_data.time_shape_id
    time_shape_id <= PULSEQ_DEFAULT_TIME_SHAPE_ID && return Δt_rf / 2
    return first(decompress_shape(shape_library[time_shape_id]...)) * Δt_rf
end

function emit_gradients_section!(io::IO, event_libraries)
    write(io, "\n# Format of arbitrary gradient events:\n")
    write(io, "# id      amp      first      last  shape_id  time_id  delay\n")
    write(io, "# ..     Hz/m       Hz/m      Hz/m        ..       ..     us\n")
    write(io, "[GRADIENTS]\n")
    for id in dense_event_ids(event_libraries.grad_library)
        grad_data = event_libraries.grad_library[id]
        grad_data isa PulseqArbGradEvent || continue
        @printf(io, "%d ", id)
        emit_pulseq_float!(io, γ * grad_data.amplitude); write(io, " ")
        emit_pulseq_float!(io, γ * grad_data.first); write(io, " ")
        emit_pulseq_float!(io, γ * grad_data.last)
        @printf(io, " %d %d %d\n", grad_data.amp_shape_id, grad_data.time_shape_id, round(Int, grad_data.delay * 1e6))
    end
end

function emit_trap_section!(io::IO, event_libraries)
    write(io, "\n# Format of trapezoid gradient events:\n")
    write(io, "# id      amp      rise  flat  fall  delay\n")
    write(io, "# ..     Hz/m        us    us    us     us\n")
    write(io, "[TRAP]\n")
    for id in dense_event_ids(event_libraries.grad_library)
        trap_data = event_libraries.grad_library[id]
        trap_data isa PulseqTrapGradEvent || continue
        @printf(io, "%2d ", id)
        emit_pulseq_float!(io, γ * trap_data.amplitude)
        @printf(io, " %3d %4d %3d %3d\n", round(Int, trap_data.rise * 1e6), round(Int, trap_data.flat * 1e6), round(Int, trap_data.fall * 1e6), round(Int, trap_data.delay * 1e6))
    end
end

function emit_adc_section!(io::IO, event_libraries)
    write(io, "\n# Format of ADC events:\n")
    write(io, "# id  num  dwell  delay  freq_ppm  phase_ppm  freq  phase  phase_id\n")
    write(io, "# ..   ..     ns     us       ppm        ppm    Hz    rad        ..\n")
    write(io, "[ADC]\n")
    isempty(event_libraries.adc_library) && return
    for id in dense_event_ids(event_libraries.adc_library)
        adc_data = event_libraries.adc_library[id]
        @printf(io, "%d %d ", id, adc_data.num)
        emit_pulseq_float!(io, adc_data.dwell * 1e9); write(io, " ")
        emit_pulseq_float!(io, adc_data.delay * 1e6); write(io, " ")
        emit_pulseq_float!(io, adc_data.freq_ppm); write(io, " ")
        emit_pulseq_float!(io, adc_data.phase_ppm); write(io, " ")
        emit_pulseq_float!(io, adc_data.freq); write(io, " ")
        emit_pulseq_float!(io, adc_data.phase)
        @printf(io, " %d\n", adc_data.phase_id)
    end
end

function emit_extension_section!(io::IO, event_libraries)
    write(io, "\n# Format of extension events:\n")
    write(io, "# id  type  ref  next_id\n")
    write(io, "# Extension list is followed by extension specifications\n")
    write(io, "[EXTENSIONS]\n")
    isempty(event_libraries.extension_instance_library) && return
    for id in dense_event_ids(event_libraries.extension_instance_library)
        ext_data = event_libraries.extension_instance_library[id]
        @printf(io, "%d %d %d %d\n", id, ext_data.type, ext_data.ref, ext_data.next_id)
    end
    write(io, "\n")
    type_ids = dense_event_ids(event_libraries.extension_type_library)
    for id in type_ids
        ext_type = event_libraries.extension_type_library[id]
        write(io, extension_type_header(ext_type))
        write(io, "extension $(string(get_symbol_from_EXT_type(ext_type))) $id\n")
        for spec_id in dense_event_ids(event_libraries.extension_spec_library[id])
            emit_extension_spec!(io, spec_id, event_libraries.extension_spec_library[id][spec_id])
        end
        id != last(type_ids) && write(io, "\n")
    end
end

function emit_shapes_section!(io::IO, event_libraries)
    write(io, "\n# Sequence shapes\n")
    write(io, "[SHAPES]\n\n")
    for id in dense_event_ids(event_libraries.shape_library)
        num_samples, samples = event_libraries.shape_library[id]
        write(io, "shape_id $id\n")
        write(io, "num_samples $num_samples\n")
        for sample in samples
            emit_pulseq_float!(io, sample)
            write(io, "\n")
        end
        write(io, "\n")
    end
end

function emit_signature_section!(io, algorithm, hash_value)
    write(io, "[SIGNATURE]\n")
    write(io, "# This is the hash of the Pulseq file, calculated right before the [SIGNATURE] section was added\n")
    write(io, "# It can be reproduced/verified with md5sum if the file trimmed to the position right above [SIGNATURE]\n")
    write(io, "# The new line character preceding [SIGNATURE] BELONGS to the signature (and needs to be stripped away for recalculating/verification)\n")
    write(io, "Type $(algorithm)\n")
    write(io, "Hash $(hash_value)\n\n")
end

const INTERNAL_DEFINITION_KEYS = ("signature", "PulseqVersion", "FileName", KomaMRIBase.PULSEQ_HW_DEFINITION_KEYS...)
const PULSEQ_RASTER_DEFINITION_KEYS = KomaMRIBase.PULSEQ_RASTER_DEFINITION_KEYS

should_emit_definition(key, value) =
    key ∉ INTERNAL_DEFINITION_KEYS &&
    key ∉ PULSEQ_RASTER_DEFINITION_KEYS &&
    key != "RequiredExtensions" &&
    !(value isa AbstractVector && isempty(value))

emit_definition_value!(io, value::AbstractString) = write(io, value)
emit_definition_value!(io, value::Char) = write(io, string(value))
emit_definition_value!(io, value::Real) = emit_pulseq_float!(io, value)

function emit_definition_value!(io, values::Union{AbstractVector,Tuple})
    first_value = true
    for value in values
        first_value || write(io, " ")
        emit_definition_value!(io, value)
        first_value = false
    end
end

emit_definition_value!(io, value) = write(io, string(value))

definition_entry_value(value::Union{AbstractVector,Tuple}) = join(definition_entry_value.(value), " ")
definition_entry_value(value::Real) = pulseq_float_string(value)
definition_entry_value(value) = string(value)

function pulseq_definitions(definitions, raster)
    entries = Pair{String,String}[]
    required_extensions = String[]
    for (key, value) in definitions
        if key == "RequiredExtensions"
            required_extensions = value isa AbstractVector ? string.(value) : [string(value)]
        elseif should_emit_definition(key, value)
            push!(entries, key => definition_entry_value(value))
        end
    end
    return PulseqDefinitions(
        entries,
        raster.BlockDurationRaster,
        raster.GradientRasterTime,
        raster.RadiofrequencyRasterTime,
        raster.AdcRasterTime,
        required_extensions,
    )
end

function emit_extension_spec!(io, spec_id, spec::Trigger)
    @printf(io, "%d %d %d %d %d\n", spec_id, spec.type, spec.channel, round(Int, spec.delay * 1e6), round(Int, spec.duration * 1e6))
end

function emit_extension_spec!(io, spec_id, spec::Union{LabelSet,LabelInc})
    @printf(io, "%d %d %s\n", spec_id, spec.labelvalue, spec.labelstring)
end

function emit_extension_spec!(io, spec_id, spec::Extension)
    write(io, string(spec_id))
    for field in fieldnames(typeof(spec))
        write(io, " ")
        emit_definition_value!(io, getfield(spec, field))
    end
    write(io, "\n")
end

function prepare_pulseq_write(seq, raster=DEFAULT_RASTER)
    return prepare_pulseq_write(seq.GR, seq.RF, seq.ADC, seq.DUR, seq.EXT, seq.DEF, raster)
end

function prepare_pulseq_write(gr, rf, adc, block_durations, ext, def, raster)
    warn_count = Threads.Atomic{Int}(0)
    rf_type = eltype(rf)
    grad_type = eltype(gr)
    canonical_rfs = IdDict{rf_type,rf_type}()
    canonical_grads = IdDict{grad_type,grad_type}()
    canonical_adcs = IdDict{ADC,ADC}()
    axis_names = ("x", "y", "z")
    qgr = Matrix{grad_type}(undef, size(gr))
    qrf = Matrix{rf_type}(undef, size(rf))
    qadc = similar(adc)
    qdur = copy(block_durations)
    n_thread_slots = Threads.maxthreadid()
    thread_rf_caches = [IdDict{rf_type,rf_type}() for _ in 1:n_thread_slots]
    thread_grad_caches = [IdDict{grad_type,grad_type}() for _ in 1:n_thread_slots]
    thread_adc_caches = [IdDict{ADC,Tuple{ADC,Float64}}() for _ in 1:n_thread_slots]
    Threads.@threads for bi in eachindex(block_durations)
        qrf[1, bi], qgr[1, bi], qgr[2, bi], qgr[3, bi], qadc[bi], qdur[bi] =
            quantize_block(
                rf,
                gr,
                adc,
                block_durations,
                raster,
                axis_names,
                bi,
                warn_count,
                thread_rf_caches[Threads.threadid()],
                thread_grad_caches[Threads.threadid()],
                thread_adc_caches[Threads.threadid()],
            )
    end
    canonicalize_quantized!(qrf, canonical_rfs)
    canonicalize_quantized!(qgr, canonical_grads)
    canonicalize_quantized!(qadc, canonical_adcs)
    if any(block -> any(ext -> ext isa QuaternionRot, block), ext)
        qgr = KomaMRIBase._apply_rotation_extensions_to_gradients!(qgr, ext; reverse=true)
    end
    return KomaMRIBase.Sequence(qgr, qrf, qadc, qdur, ext, def)
end

function quantize_block(
    rf,
    gr,
    adc,
    block_durations,
    raster,
    axis_names,
    bi,
    warn_count,
    rf_cache,
    grad_cache,
    adc_cache,
)
    rf_i = rf[1, bi]
    qrf = get_quantized!(rf_cache, rf_i) do
        is_on(rf_i) || return rf_i
        qrfi = copy(rf_i)
        quantize_rf!(qrfi, raster, bi, warn_count)
    end
    qgr = ntuple(Val(3)) do gi
        grad_i = gr[gi, bi]
        get_quantized!(grad_cache, grad_i) do
            is_on(grad_i) || return grad_i
            qgrad_i = copy(grad_i)
            quantize_grad!(qgrad_i, raster, bi, axis_names[gi], warn_count)
        end
    end

    adc_i = adc[bi]
    qadc, adc_end = get_quantized!(adc_cache, adc_i) do
        is_ADC_on(adc_i) || return adc_i, 0.0
        qadci = copy(adc_i)
        qadci, quantize_adc!(qadci, raster, bi, warn_count)
    end

    key = :BlockDurationRaster
    max_event_duration = max(dur(qgr[1]), dur(qgr[2]), dur(qgr[3]), dur(qrf), adc_end)
    block_duration = max(block_durations[bi], max_event_duration)
    warn_ctx = (; block_id=bi, event=:Block, warn_count)
    qdur = quantize_time_value(block_duration, getfield(raster, key), warn_ctx, :DUR; step_name=key)
    return qrf, qgr..., qadc, qdur
end

function quantize_rf!(rf, raster, block_id, warn_count)
    key = :RadiofrequencyRasterTime
    Δt_rf = getfield(raster, key)
    warn_ctx = (; block_id, event=:RF, warn_count)
    qt(value, field; step=Δt_rf, step_name=key) =
        quantize_time_value(value, step, warn_ctx, field; step_name)
    qT(value, field; step=Δt_rf, step_name=key) =
        quantize_sample_times(value, sample_count(rf), step, warn_ctx, field; step_name)
    policy = pulseq_rf_timing_policy(rf, Δt_rf)
    rf.T = qT(rf.T, :T; step=policy.time_step)
    policy = pulseq_rf_timing_policy(rf, Δt_rf)
    if isnothing(policy.time_shape_id)
        rf.delay = qt(rf.delay, :delay)
    else
        pulseq_delay = max(rf.delay - policy.first_sample_offset, 0.0)
        pulseq_delay = qt(pulseq_delay, :pulseq_delay)
        rf.delay = pulseq_delay + policy.first_sample_offset
    end
    return rf
end

uniform_pulseq_interval(::Grad) = nothing

function uniform_pulseq_interval(gr::UniformlySampledGrad)
    n = length(gr.A)
    return n > 1 ? gr.T / (n - 1) : nothing
end

function uniform_pulseq_interval(gr::TimeShapedGrad)
    isempty(gr.T) && return nothing
    interval = gr.T[1]
    all(isapprox(Δt, interval; rtol=0, atol=PULSEQ_TIME_TOL) for Δt in gr.T) || return nothing
    return interval
end

function pulseq_grad_compact_time_shape_id(gr, Δt_gr)
    interval = uniform_pulseq_interval(gr)
    isnothing(interval) && return nothing
    isapprox(gr.rise, Δt_gr / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    isapprox(gr.fall, Δt_gr / 2; rtol=0, atol=PULSEQ_TIME_TOL) || return nothing
    if isapprox(interval, Δt_gr; rtol=0, atol=PULSEQ_TIME_TOL)
        return PULSEQ_DEFAULT_TIME_SHAPE_ID
    elseif isodd(length(gr.A)) && isapprox(interval, Δt_gr / 2; rtol=0, atol=PULSEQ_TIME_TOL)
        return PULSEQ_OVERSAMPLED_TIME_SHAPE_ID
    end
    return nothing
end

function pulseq_grad_timing_policy(::Grad, Δt_gr)
    return (; time_shape_id=nothing, time_step=Δt_gr, rise_fall_step=Δt_gr)
end

function pulseq_grad_timing_policy(gr::Union{UniformlySampledGrad,TimeShapedGrad}, Δt_gr)
    time_shape_id = pulseq_grad_compact_time_shape_id(gr, Δt_gr)
    time_step = time_shape_id == PULSEQ_OVERSAMPLED_TIME_SHAPE_ID ? Δt_gr / 2 : Δt_gr
    rise_fall_step = isnothing(time_shape_id) ? Δt_gr : Δt_gr / 2
    return (; time_shape_id, time_step, rise_fall_step)
end

function quantize_grad!(gr, raster, block_id, axis_name, warn_count)
    key = :GradientRasterTime
    event = Symbol("GR", axis_name)
    Δt_gr = getfield(raster, key)
    warn_ctx = (; block_id, event, warn_count)
    qt(value, field; step=Δt_gr, step_name=key) =
        quantize_time_value(value, step, warn_ctx, field; step_name)
    qT(value, field; step=Δt_gr, step_name=key) =
        quantize_sample_times(value, sample_count(gr), step, warn_ctx, field; step_name)
    policy = pulseq_grad_timing_policy(gr, Δt_gr)
    gr.delay = qt(gr.delay, :delay)
    gr.T     = qT(gr.T, :T; step=policy.time_step)
    gr.rise  = qt(gr.rise, :rise; step=policy.rise_fall_step)
    gr.fall  = qt(gr.fall, :fall; step=policy.rise_fall_step)
    return gr
end

function quantize_adc!(adc, raster, block_id, warn_count)
    warn_ctx = (; block_id, event=:ADC, warn_count)
    qt(value, field; step, step_name) =
        quantize_time_value(value, step, warn_ctx, field; step_name)
    dwell = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
    dwell = qt(dwell, :dwell_time; step=raster.AdcRasterTime, step_name=:AdcRasterTime)
    adc.T = adc.N == 1 ? dwell : (adc.N - 1) * dwell
    pulseq_delay = adc.delay - dwell / 2
    if pulseq_delay < 0
        @warn "ADC delay in Koma ($(round(adc.delay*1e3, digits=4)) ms) is below dwell/2 ($(round(dwell/2*1e3, digits=4)) ms).
In Pulseq, the first ADC sample is acquired at dwell/2. Therefore, a Pulseq delay of 0 corresponds to a Koma delay of dwell/2.
This means the Koma delay must be >= dwell/2. Clamping it to this minimum value...
See https://pulseq.github.io/pulseq_shapes_and_times.pdf#page=10"
        pulseq_delay = 0
    end
    pulseq_delay = qt(pulseq_delay, :pulseq_delay; step=raster.RadiofrequencyRasterTime, step_name=:RadiofrequencyRasterTime)
    adc.delay = pulseq_delay + dwell / 2
    return adc.delay + adc.T + dwell / 2
end

function quantize_time_value(t::Number, step, warn_ctx, field; step_name=nothing)
    rounded = round(t / step) * step
    if isapprox(t, rounded; rtol=0, atol=PULSEQ_TIME_TOL)
        return rounded
    end
    warn_time_quantization(warn_ctx, field, step; step_name)
    return ceil(t / step - PULSEQ_TIME_TOL / step) * step
end

sample_count(::Union{BlockPulseRF,TrapezoidalGrad}) = 2
sample_count(x::Union{RF,Grad}) = length(x.A)

function quantize_time_value(t::AbstractVector{<:Real}, step, warn_ctx, field; step_name=nothing)
    rounded = round.(t ./ step) .* step
    if all(isapprox(t[i], rounded[i]; rtol=0, atol=PULSEQ_TIME_TOL) for i in eachindex(t, rounded))
        return rounded
    end
    warn_time_quantization(warn_ctx, field, step; step_name)
    return ceil.(t ./ step .- PULSEQ_TIME_TOL / step) .* step
end

quantize_sample_times(t::AbstractVector{<:Real}, n_samples, step, warn_ctx, field; step_name=nothing) =
    quantize_time_value(t, step, warn_ctx, field; step_name)

function quantize_sample_times(t::Number, n_samples, step, warn_ctx, field; step_name=nothing)
    n_intervals = max(n_samples - 1, 1)
    interval = n_intervals == 1 ? t : t / n_intervals
    q_interval = quantize_time_value(interval, step, warn_ctx, field; step_name)
    return n_intervals == 1 ? q_interval : q_interval * n_intervals
end

warning_label(x::Symbol) = replace(String(x), '_' => ' ')
warning_label(x) = string(x)

pretty_time_step(step) = step < 1e-6 ? "$(step * 1e9) ns" : "$(step * 1e6) μs"

function warn_time_quantization(warn_ctx, field, step; step_name=nothing)
    count = Threads.atomic_add!(warn_ctx.warn_count, 1)
    if count < 10
        @warn "Block $(warn_ctx.block_id): Event $(warning_label(warn_ctx.event)): Element $(warning_label(field)):
Time is not a multiple of $(isnothing(step_name) ? pretty_time_step(step) : "$(warning_label(step_name)) ($(pretty_time_step(step)))"). Quantizing it..."
    elseif count == 10
        @warn "Additional time quantization warnings occurred; detailed logs are capped at 10."
    end
    return nothing
end

"""
    write_seq(seq, filename)

Writes a Sequence struct to a Pulseq file with `.seq` extension.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `filename`: (`::String`) absolute or relative path of the sequence file `.seq`
- `sys`: optional scanner used for hardware-limit validation and export raster
  times. If passed, `sys` raster times take precedence over `seq.DEF`.
- `verbose`: (`::Bool`, `=true`) show informational writing/checking messages

# Examples
```julia-repl
julia> seq = Sequence()
julia> @addblock seq += RF(2e-6, 3e-3)
julia> @addblock seq += ADC(100, 100e-3)
julia> write_seq(seq, "fid.seq")
```
"""
function write_seq(
    data::PulseqSequenceData, filename::AbstractString;
    signatureAlgorithm="md5",
    verbose=true,
)
    verbose && @info "Saving sequence to $(basename(filename)) ..."
    payload = let io = IOBuffer()
        emit_pulseq(io, data)
        take!(io)
    end
    signature_hash = supported_signature_digest(signatureAlgorithm, payload)
    open(filename, "w") do io
        write(io, payload)
        write(io, '\n')
        emit_signature_section!(io, signatureAlgorithm, signature_hash)
    end
    return nothing
end

function pulseq_data(seq::KomaMRIBase.Sequence; sys=nothing, check_timing=true, check_hw_limits=true, verbose=true)
    if check_hw_limits
        isnothing(sys) ? KomaMRIBase.check_hw_limits(seq) : KomaMRIBase.check_hw_limits(seq, sys)
    end
    raster = isnothing(sys) ? PulseqRaster(seq) : PulseqRaster(seq, sys)
    if verbose
        @info "Checking Pulseq raster..." BlockDurationRaster = raster.BlockDurationRaster GradientRasterTime = raster.GradientRasterTime RadiofrequencyRasterTime = raster.RadiofrequencyRasterTime AdcRasterTime = raster.AdcRasterTime
    end
    prepared = prepare_pulseq_write(seq, raster)
    if check_timing
        KomaMRIBase.check_timing(prepared, pulseq_raster_tuple(raster))
    end
    blocks, event_libraries = collect_pulseq_assets(prepared, raster)
    return PulseqSequenceData(blocks, event_libraries, v"1.5.1", nothing)
end

"""
    data = write_seq_data(seq; sys=nothing, check_timing=true, check_hw_limits=true, verbose=true)

Prepare a Koma `Sequence` as [`PulseqSequenceData`](@ref) without writing a file.
Without `sys`, timing and hardware checks use metadata stored in `seq.DEF`. With
`sys`, timing, hardware checks, and export raster times use `sys`.
"""
write_seq_data(seq::KomaMRIBase.Sequence; sys=nothing, check_timing=true, check_hw_limits=true, verbose=true) =
    pulseq_data(seq; sys, check_timing, check_hw_limits, verbose)

function write_seq(
    seq::KomaMRIBase.Sequence, filename::AbstractString;
    sys=nothing,
    signatureAlgorithm="md5",
    check_timing=true,
    check_hw_limits=true,
    verbose=true,
)
    data = pulseq_data(seq; sys, check_timing, check_hw_limits, verbose)
    return write_seq(data, filename; signatureAlgorithm, verbose)
end

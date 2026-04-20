"""
    collect_pulseq_assets(seq::Sequence, raster::PulseqRaster) -> (blocks, event_libraries)

Create the Pulseq export dictionaries required to serialize `seq` into the Pulseq file format.
This function is responsible for deduplicating reusable objects (RF, gradients, shapes, etc.)
and for translating each sequence block into integer lookups expected by the specification.
"""
function collect_pulseq_assets(seq::Sequence, raster::PulseqRaster)
    blocks = NTuple{7,Int}[]
    rf_library = Dict{Int,RFEvent}()
    grad_library = Dict{Int,GradEvent}()
    adc_library = Dict{Int,ADCEvent}()
    tmp_delay_library = Dict{Int,Float64}()
    shape_library = ShapeLibrary()
    extension_instance_library = Dict{Int,ExtensionInstanceEvent}()
    extension_type_library = Dict{Int,Type{<:Extension}}()
    extension_spec_library = Dict{Int,Dict{Int,Extension}}()
    definitions = Dict{String, Any}(seq.DEF)

    # Hash indexes to avoid linear scans while deduplicating events/shapes.
    rf_index = Dict{UInt, Vector{Int}}()
    grad_index = Dict{UInt, Vector{Int}}()
    adc_index = Dict{UInt, Vector{Int}}()
    ext_instance_index = Dict{UInt, Vector{Int}}()
    ext_type_index = Dict{UInt, Vector{Int}}()
    ext_spec_index = Dict{Int, Dict{UInt, Vector{Int}}}()
    shape_exact_index = Dict{UInt, Vector{Int}}()
    shape_phase_index = Dict{Int, Vector{Int}}()
    
    # First step: collect all unique extension vectors and register them once.
    extension_vectors = Vector{Extension}[]
    extension_vector_ids = Int[]

    for block in seq
        ext_vec = block.EXT[1]
        isempty(ext_vec) && continue
        matching_idx = findfirst(
            registered_ext ->
                length(registered_ext) == length(ext_vec) &&
                all(e1 ≈ e2 for (e1, e2) in zip(registered_ext, ext_vec)),
            extension_vectors,
        )
        if matching_idx === nothing
            ext_id = register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext_vec; 
                ext_instance_index=ext_instance_index, ext_type_index=ext_type_index, ext_spec_index=ext_spec_index)
            push!(extension_vectors, ext_vec)
            push!(extension_vector_ids, ext_id)
        end
    end

    # Second step: process all blocks and use pre-registered extension IDs.
    for block in seq
        rf = block.RF[1]
        rf_id = register_rf!(rf_library, shape_library, rf.A, rf.T, rf.Δf, rf.delay, rf.center, rf.use, raster;
            rf_index=rf_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)

        gx, gy, gz = block.GR
        gx_id = register_grad!(
            grad_library, shape_library, gx.A, gx.T, gx.rise, gx.fall, gx.delay, gx.first, gx.last, raster;
            grad_index=grad_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
        )
        gy_id = register_grad!(
            grad_library, shape_library, gy.A, gy.T, gy.rise, gy.fall, gy.delay, gy.first, gy.last, raster;
            grad_index=grad_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
        )
        gz_id = register_grad!(
            grad_library, shape_library, gz.A, gz.T, gz.rise, gz.fall, gz.delay, gz.first, gz.last, raster;
            grad_index=grad_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
        )

        adc = block.ADC[1]
        adc_id, adc_dur = register_adc!(adc_library, adc.N, adc.T, adc.delay, adc.Δf, adc.ϕ; adc_index=adc_index)

        max_event_duration = max(dur(gx), dur(gy), dur(gz), dur(rf), adc_dur)
        block_duration = max(dur(block), max_event_duration)
        ratio = block_duration / raster.BlockDurationRaster
        base = floor(Int, ratio)
        frac = ratio - base
        round_up_threshold = 0.1
        duration = base + (frac >= round_up_threshold - QUANT_TOL ? 1 : 0)

        ext_vec = block.EXT[1]
        ext_id = 0
        if !isempty(ext_vec)
            matching_idx = findfirst(
                registered_ext ->
                    length(registered_ext) == length(ext_vec) &&
                    all(e1 ≈ e2 for (e1, e2) in zip(registered_ext, ext_vec)),
                extension_vectors,
            )
            ext_id = matching_idx === nothing ? 0 : extension_vector_ids[matching_idx]
        end
        push!(blocks, (duration, rf_id, gx_id, gy_id, gz_id, adc_id, ext_id))
    end

    merge_definitions_with_raster!(definitions, raster)

    event_libraries = PulseqEventLibraries(
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
    id = register_rf!(rf_library, shape_library, A, T, Δf, delay, center, use, raster)
"""
function register_rf!(
    rf_library, shape_library, A, T, Δf, delay, center, use, raster::PulseqRaster;
    rf_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    iszero(maximum(abs.(A))) && return 0
    Δt_rf = raster.RadiofrequencyRasterTime
    mag_id, phase_id, time_id = register_rf_shapes!(shape_library, A, T, Δt_rf;
        compress=true, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
    )
    amp = γ * maximum(abs.(A)) # from T to Hz (nucleus-dependent)
    delay = delay - (time_id == 0) * Δt_rf / 2
    phase  = mod.(angle.(A), 2π)[1]
    use_char = KomaMRIBase.get_char_from_RF_use(use)
    rf_event = RFEvent(amp, mag_id, phase_id, time_id, center, delay, 0.0, 0.0, Δf, phase, use_char)
    return _store_event!(rf_library, rf_event, rf_index)
end

# A and T are numbers (pulse waveform)
function register_rf_shapes!(
    shapes, A, T, Δrf; compress = true, 
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    mag_id   = _store_shape!(shapes, [1.0, 1.0]; compress=compress, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    phase_id = _store_shape!(shapes, [0.0, 0.0]; compress=compress, phase_shape=true, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    time_id  = _store_shape!(shapes, [0.0, T] ./ Δrf; compress=compress, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    return mag_id, phase_id, time_id
end
# A is a vector (uniformly-sampled waveform)
function register_rf_shapes!(
    shapes, A::Vector, T, Δrf; compress = true,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    n_samples   = length(A)
    mag_id      = _store_shape!(shapes, abs.(A) ./ maximum(abs.(A)); compress=compress, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    phase_shape = mod.(angle.(A), 2π) / 2π
    phase_id    = _store_shape!(shapes, phase_shape .- phase_shape[1]; compress=compress, phase_shape=true, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    t_vector    = collect(range(0, T, length=n_samples)) ./ Δrf
    time_id     = _store_shape!(shapes, t_vector; compress=compress, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    return mag_id, phase_id, time_id
end
# A and T are vectors (time-shaped waveform)
function register_rf_shapes!(
    shapes, A::Vector, T::Vector, Δrf; compress = true,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    mag_id      = _store_shape!(shapes, abs.(A) ./ maximum(abs.(A)); compress=compress, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    phase_shape = mod.(angle.(A), 2π) / 2π
    phase_id    = _store_shape!(shapes, phase_shape .- phase_shape[1]; compress=compress, phase_shape=true, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    t_vector    = cumsum(vcat(0.0, T ./ Δrf))
    time_id     = _store_shape!(shapes, t_vector; compress=compress, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
    return mag_id, phase_id, time_id
end

"""
    id = register_grad!(grad_library, shape_library, A, T, rise, fall, delay, first, last, raster)
"""
# Time-shaped waveform
function register_grad!(
    grad_library, shape_library, A::Vector, T::Vector, rise, fall, delay, first, last, raster::PulseqRaster;
    grad_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    iszero(maximum(abs.(A))) && return 0
    if (iszero(rise) && iszero(fall))
        shape_id = _store_shape!(shape_library, A ./ maximum(abs.(A)); compress=true, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
        t_vector = cumsum([0; T]) ./ raster.GradientRasterTime
        time_id  = _store_shape!(shape_library, t_vector; compress=true, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index)
        amp   = γ * maximum(abs.(A)) # from T/m to Hz/m (nucleus-dependent)
        first = γ * first # from T/m to Hz/m 
        last  = γ * last  # from T/m to Hz/m
        return _store_event!(grad_library, ArbGradEvent(amp, first, last, shape_id, time_id, delay), grad_index)
    end
    return register_grad!(
        grad_library, shape_library, [first; A; last], [rise; T; fall], 0, 0, delay, first, last, raster;
        grad_index=grad_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
    )
end
# Uniformly-sampled waveform
function register_grad!(
    grad_library, shape_library, A::Vector, T::Number, rise, fall, delay, first, last, raster::PulseqRaster;
    grad_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    iszero(maximum(abs.(A))) && return 0
    intervals = diff(collect(range(0, T, length=length(A))))
    return register_grad!(
        grad_library, shape_library, [first; A; last], [rise; intervals; fall], 0, 0, delay, first, last, raster;
        grad_index=grad_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
    )
end
# Trapezoidal waveform
function register_grad!(
    grad_library, shape_library, A::Number, T::Number, rise, fall, delay,  first, last,  raster::PulseqRaster;
    grad_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    iszero(A) && return 0
    if (iszero(first) && iszero(last)) 
        return _store_event!(grad_library, TrapGradEvent(γ * A, rise, T, fall, delay), grad_index)
    end
    return register_grad!(
        grad_library, shape_library, [first, A, A, last], [rise, T, fall], 0, 0, delay, first, last, raster;
        grad_index=grad_index, shape_exact_index=shape_exact_index, shape_phase_index=shape_phase_index
    )
end 

"""
    id = register_adc!(adc_library, N, T, delay, Δf, ϕ)

In Pulseq the ADC event starts at a dwell-time edge, but the first sample is taken at dwell/2
(see [Pulseq time and shape specification](https://pulseq.github.io/pulseq_shapes_and_times.pdf)).
Koma uses delay = time to first sample, so we store pulseq_delay = delay - dwell_s/2. For exact
round-trip (write → read) the Koma ADC delay must be ≥ dwell_s/2; otherwise it is clamped and a
warning is emitted.
"""
function register_adc!(adc_library, N, T, delay, Δf, ϕ; adc_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing)
    iszero(N) && return (0, 0.0)
    dwell_s = N == 1 ? T : T / (N - 1)
    delay_s = delay - dwell_s / 2
    adc_event = ADCEvent(N, dwell_s, delay_s, 0.0, 0.0, Δf, ϕ, 0)
    return _store_event!(adc_library, adc_event, adc_index), delay_s + T + dwell_s
end

"""
    id = register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext)
"""
function register_ext!(
    extension_instance_library, extension_type_library, extension_spec_library, ext::Vector{Extension};
    ext_instance_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    ext_type_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    ext_spec_index::Union{Nothing, Dict{Int, Dict{UInt, Vector{Int}}}}=nothing,
)
    length(ext) == 0 && return 0
    instance_ids = Int[]
    for e in ext
        ext_id = _store_event!(extension_type_library, typeof(e), ext_type_index)
        if !haskey(extension_spec_library, ext_id)
            extension_spec_library[ext_id] = Dict{Int,Extension}()
            !isnothing(ext_spec_index) && (ext_spec_index[ext_id] = Dict{UInt, Vector{Int}}())
        end
        ref = _store_event!(
            extension_spec_library[ext_id],
            e,
            isnothing(ext_spec_index) ? nothing : get!(ext_spec_index, ext_id, Dict{UInt, Vector{Int}}()),
        )
        instance_id = _store_event!(extension_instance_library, ExtensionInstanceEvent((ext_id, ref, 0)), ext_instance_index)
        push!(instance_ids, instance_id)
    end

    for i in 1:(length(instance_ids) - 1)
        current_id = instance_ids[i]
        next_id = instance_ids[i + 1]
        current = extension_instance_library[current_id]
        extension_instance_library[current_id] = ExtensionInstanceEvent(current.type, current.ref, next_id)
    end

    return instance_ids[1]
end

# store_event! and store_shape! are helper functions to store events and shapes in the dictionary, and return the id of the event or shape.
function _store_event!(event_dict::Dict, event)
    for (k, v) in event_dict
        if v == event return k end
    end 
    new_id = maximum(keys(event_dict); init=0) + 1
    event_dict[new_id] = event
    return new_id
end

function _store_event!(event_dict::Dict, event, event_index::Nothing)
    return _store_event!(event_dict, event)
end

function _store_event!(event_dict::Dict, event, event_index::Dict{UInt, Vector{Int}})
    h = hash(event)
    for id in get(event_index, h, Int[])
        if event_dict[id] == event
            return id
        end
    end
    new_id = length(event_dict) + 1
    event_dict[new_id] = event
    push!(get!(event_index, h, Int[]), new_id)
    return new_id
end

@inline _shape_hash(num_samples::Int, samples::Vector{Float64}) = hash(samples, hash(num_samples))

function _store_shape!(
    shapes::Dict{Int,Tuple{Int,Vector{Float64}}},
    samples::Vector{Float64};
    compress::Bool = true,
    phase_shape::Bool = false,
    shape_exact_index::Union{Nothing, Dict{UInt, Vector{Int}}}=nothing,
    shape_phase_index::Union{Nothing, Dict{Int, Vector{Int}}}=nothing,
)
    payload = compress ? compress_shape(samples) : (length(samples), samples)
    isnothing(shape_exact_index) && isnothing(shape_phase_index) && begin
        for (k, existing) in shapes
            if phase_shape
                existing_phase = existing[2] .- existing[2][1]
                if safe_equal_angles(existing_phase, payload[2]) && existing[1] == payload[1] return k end
            else
                if existing == payload return k end
            end
        end
        new_id = length(shapes) + 1
        shapes[new_id] = payload
        return new_id
    end
    if phase_shape
        for id in get(shape_phase_index, payload[1], Int[])
            existing = shapes[id]
            existing_phase = existing[2] .- existing[2][1]
            if safe_equal_angles(existing_phase, payload[2]) && existing[1] == payload[1]
                return id
            end
        end
        new_id = length(shapes) + 1
        shapes[new_id] = payload
        push!(get!(shape_phase_index, payload[1], Int[]), new_id)
        return new_id
    end
    h = _shape_hash(payload[1], payload[2])
    for id in get(shape_exact_index, h, Int[])
        if shapes[id] == payload
            return id
        end
    end
    new_id = length(shapes) + 1
    shapes[new_id] = payload
    push!(get!(shape_exact_index, h, Int[]), new_id)
    return new_id
end

safe_equal_angles(a1, a2) = length(a1) == length(a2) && abs(sum(exp.(1im * 2π * a1) .* exp.(-1im * 2π * a2)) / length(a1)) == 1

"""
    emit_pulseq(io::IO, blocks, event_libraries)

Write the Pulseq sections to `io`, using the already prepared `blocks` and `event_libraries`.
"""
function emit_pulseq(io::IO, blocks, event_libraries)
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
    if any(g -> g isa ArbGradEvent, values(event_libraries.grad_library))
        emit_gradients_section!(io, event_libraries)
    end
    if any(g -> g isa TrapGradEvent, values(event_libraries.grad_library))
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

emit_header_comment!(io::IO) = write(io, "# Pulseq sequence file\n# Created by KomaMRI.jl\n\n")

function emit_version_section!(io::IO, version=v"1.5.1")
    write(io, "[VERSION]\n")
    write(io, "major $(version.major)\n")
    write(io, "minor $(version.minor)\n")
    write(io, "revision $(version.patch)\n\n")
end

function emit_definitions_section!(io::IO, definitions)
    write(io, "[DEFINITIONS]\n")
    for (key, value) in definitions
        write(io, "$key $value\n")
    end
    write(io, "\n")
end

# --- Pulseq table formatting ---
const PULSEQ_TABLE_COL_GAP = " "

_format_value(val::AbstractFloat) =
    isapprox(val, round(val); atol=QUANT_TOL, rtol=sqrt(eps(Float64))) ? string(round(Int, val)) : @sprintf("%.6f", val)
_format_value(val) = string(val)

function _compute_col_widths(rows; header_rows=())
    ncols = !isempty(rows) ? length(first(rows)) : !isempty(header_rows) ? length(first(header_rows)) : 0
    widths = zeros(Int, ncols)
    for header in header_rows
        for i in 1:ncols
            widths[i] = max(widths[i], length(string(header[i])))
        end
    end
    for row in rows
        for i in 1:ncols
            widths[i] = max(widths[i], length(_format_value(row[i])))
        end
    end
    return widths
end

"""Reserve two extra characters in column 1 so `# ` + cell (width `w[1]-2`) matches data column width `w[1]`."""
function _pulseq_reserve_hash_prefix_width!(widths::AbstractVector{Int})
    isempty(widths) && return widths
    widths[1] += 2
    return widths
end

function _format_row(io::IO, values, max_lengths::AbstractVector{<:Integer})
    for (i, (val, max_len)) in enumerate(zip(values, max_lengths))
        str = _format_value(val)
        num_spaces = max(0, max_len - length(str))
        write(io, repeat(" ", num_spaces))
        write(io, str)
        if i < length(values)
            write(io, PULSEQ_TABLE_COL_GAP)
        end
    end
    write(io, "\n")
end

function _format_comment_row(io::IO, values, max_lengths::AbstractVector{<:Integer})
    lengths = collect(max_lengths)
    if !isempty(lengths)
        lengths[1] = max(1, lengths[1] - 2)
    end
    write(io, "# ")
    _format_row(io, values, lengths)
end

"""One extension specification line (Pulseq file units), using the same layout as read: `"%i " * get_scanf_format(T)`."""
function _pulseq_extension_spec_line(ext_type::Type{<:Extension}, spec_id::Int, v::Extension)
    fmt = strip("%i " * strip(KomaMRIBase.get_scanf_format(ext_type)))
    tokens = split(fmt, r"\s+"; keepempty=false)
    field_values = Tuple(getfield(v, f) for f in fieldnames(typeof(v)))
    scales = vec(collect(KomaMRIBase.get_scale(typeof(v))))
    vals = Any[spec_id]
    for (s, d) in zip(scales, field_values)
        push!(vals, d isa Number ? d / s : d)
    end
    length(tokens) == length(vals) ||
        error("Pulseq export: extension $ext_type: $(length(tokens)) scanf tokens but $(length(vals)) values; check `get_scanf_format` vs struct fields.")
    parts = String[]
    for (tok, x) in zip(tokens, vals)
        t = strip(tok)
        cell = if startswith(t, "%i") || startswith(t, "%d")
            @sprintf("%d", round(Int, x isa Integer ? x : Float64(x)))
        elseif startswith(t, "%f") || startswith(t, "%e") || startswith(t, "%g")
            @sprintf("%.6f", Float64(x))
        elseif startswith(t, "%s")
            string(x)
        elseif startswith(t, "%c")
            string(x isa Char ? x : first(string(x)))
        else
            string(x)
        end
        push!(parts, cell)
    end
    return join(parts, ' ') * "\n"
end

function emit_blocks_section!(io::IO, blocks)
    write(io, "# Format of blocks:\n")
    rows = NTuple{8,Any}[]
    for (id, block) in enumerate(blocks)
        duration, rf, gx, gy, gz, adc, ext = block
        push!(rows, (id, duration, rf, gx, gy, gz, adc, ext))
    end
    header_rows = (("NUM", "DUR", "RF", "GX", "GY", "GZ", "ADC", "EXT"),)
    widths = _compute_col_widths(rows; header_rows=header_rows)
    _pulseq_reserve_hash_prefix_width!(widths)
    _format_comment_row(io, header_rows[1], widths)
    write(io, "[BLOCKS]\n")
    for row in rows
        _format_row(io, row, widths)
    end
end

function emit_rf_section!(io::IO, event_libraries)
    write(io, "\n# Format of RF events:\n")
    rows = NTuple{12,Any}[]
    rf_ids = sort(collect(keys(event_libraries.rf_library)))
    for id in rf_ids
        rf_data = event_libraries.rf_library[id]
        center_us = isnothing(rf_data.center) ? 0.0 : rf_data.center * 1e6
        delay_us = round(Int, rf_data.delay * 1e6)
        push!(rows, (id, rf_data.amplitude, rf_data.mag_id, rf_data.phase_id, rf_data.time_shape_id, center_us, delay_us, rf_data.freq_ppm, rf_data.phase_ppm, rf_data.freq, rf_data.phase, rf_data.use))
    end
    header_rows = (
        ("id", "amp", "mag_id", "phase_id", "time_id", "center", "delay", "freq_ppm", "phase_ppm", "freq_off", "phase_off", "use"),
        ("..", "Hz", "..", "..", "..", "us", "us", "ppm", "rad/MHz", "Hz", "rad", "..")
    )
    widths = _compute_col_widths(rows; header_rows=header_rows)
    _pulseq_reserve_hash_prefix_width!(widths)
    _format_comment_row(io, header_rows[1], widths)
    _format_comment_row(io, header_rows[2], widths)
    write(io, "# Field 'use' is the initial of:\n")
    write(io, "# excitation refocusing inversion saturation preparation other undefined\n")
    write(io, "[RF]\n")
    for row in rows
        _format_row(io, row, widths)
    end
end

function emit_gradients_section!(io::IO, event_libraries)
    write(io, "\n# Format of arbitrary gradient events:\n")
    rows = NTuple{7,Any}[]
    grad_ids = sort(collect(keys(event_libraries.grad_library)))
    for id in grad_ids
        grad_data = event_libraries.grad_library[id]
        if grad_data isa ArbGradEvent
            push!(rows, (id, grad_data.amplitude, grad_data.first, grad_data.last, grad_data.amp_shape_id, grad_data.time_shape_id, round(Int, grad_data.delay * 1e6)))
        end
    end
    header_rows = (
        ("id", "amp", "first", "last", "shape_id", "time_id", "delay"),
        ("..", "Hz/m", "Hz/m", "Hz/m", "..", "..", "us")
    )
    widths = _compute_col_widths(rows; header_rows=header_rows)
    _pulseq_reserve_hash_prefix_width!(widths)
    _format_comment_row(io, header_rows[1], widths)
    _format_comment_row(io, header_rows[2], widths)
    write(io, "[GRADIENTS]\n")
    for row in rows
        _format_row(io, row, widths)
    end
end

function emit_trap_section!(io::IO, event_libraries)
    write(io, "\n# Format of trapezoid gradient events:\n")
    rows = NTuple{6,Any}[]
    trap_ids = sort(collect(keys(event_libraries.grad_library)))
    for id in trap_ids
        trap_data = event_libraries.grad_library[id]
        if trap_data isa TrapGradEvent
            push!(rows, (id, trap_data.amplitude, round(Int, trap_data.rise * 1e6), round(Int, trap_data.flat * 1e6), round(Int, trap_data.fall * 1e6), round(Int, trap_data.delay * 1e6)))
        end
    end
    header_rows = (
        ("id", "amp", "rise", "flat", "fall", "delay"),
        ("..", "Hz/m", "us", "us", "us", "us")
    )
    widths = _compute_col_widths(rows; header_rows=header_rows)
    _pulseq_reserve_hash_prefix_width!(widths)
    _format_comment_row(io, header_rows[1], widths)
    _format_comment_row(io, header_rows[2], widths)
    write(io, "[TRAP]\n")
    for row in rows
        _format_row(io, row, widths)
    end
end

function emit_adc_section!(io::IO, event_libraries)
    write(io, "\n# Format of ADC events:\n")
    rows = NTuple{9,Any}[]
    adc_ids = sort(collect(keys(event_libraries.adc_library)))
    for id in adc_ids
        adc_data = event_libraries.adc_library[id]
        push!(rows, (id, adc_data.num, round(Int, adc_data.dwell * 1e9), round(Int, adc_data.delay * 1e6), adc_data.freq_ppm, adc_data.phase_ppm, adc_data.freq, adc_data.phase, adc_data.phase_id))
    end
    header_rows = (
        ("id", "num", "dwell", "delay", "freq_ppm", "phase_ppm", "freq", "phase", "phase_id"),
        ("..", "..", "ns", "us", "ppm", "ppm", "Hz", "rad", "..")
    )
    widths = _compute_col_widths(rows; header_rows=header_rows)
    _pulseq_reserve_hash_prefix_width!(widths)
    _format_comment_row(io, header_rows[1], widths)
    _format_comment_row(io, header_rows[2], widths)
    write(io, "[ADC]\n")
    for row in rows
        _format_row(io, row, widths)
    end
end

function emit_extension_section!(io::IO, event_libraries)
    write(io, "\n# Format of extension events:\n")
    ext_rows = NTuple{4,Any}[]
    ext_ids = sort(collect(keys(event_libraries.extension_instance_library)))
    for id in ext_ids
        ext_data = event_libraries.extension_instance_library[id]
        push!(ext_rows, (id, ext_data.type, ext_data.ref, ext_data.next_id))
    end
    header_rows = (("id", "type", "ref", "next_id"),)
    ext_widths = _compute_col_widths(ext_rows; header_rows=header_rows)
    _pulseq_reserve_hash_prefix_width!(ext_widths)
    _format_comment_row(io, header_rows[1], ext_widths)
    write(io, "# Extension list is followed by extension specifications\n")
    write(io, "[EXTENSIONS]\n")
    isempty(ext_rows) && return
    for row in ext_rows
        _format_row(io, row, ext_widths)
    end
    write(io, "\n")
    type_ids = sort(collect(keys(event_libraries.extension_type_library)))
    for id in type_ids
        ext_type = event_libraries.extension_type_library[id]
        write(io, KomaMRIBase.extension_type_header(ext_type))
        write(io, "extension $(string(KomaMRIBase.get_symbol_from_EXT_type(ext_type))) $id\n")
        for sid in sort(collect(keys(event_libraries.extension_spec_library[id])))
            v = event_libraries.extension_spec_library[id][sid]
            write(io, _pulseq_extension_spec_line(ext_type, sid, v))
        end
        id !== type_ids[end] && write(io, "\n")
    end
end

function emit_shapes_section!(io::IO, event_libraries)
    write(io, "\n# Sequence shapes\n")
    write(io, "[SHAPES]\n\n")
    for id in sort(collect(keys(event_libraries.shape_library)))
        num_samples, samples = event_libraries.shape_library[id]
        write(io, "shape_id $id\n")
        write(io, "num_samples $num_samples\n")
        for sample in samples
            @printf(io, "%.8f\n", sample)
        end
        write(io, "\n")
    end
end

function emit_signature_section!(io::IO, algorithm::AbstractString, hash_value::AbstractString)
    write(io, "[SIGNATURE]\n")
    write(io, "# This is the hash of the Pulseq file, calculated right before the [SIGNATURE] section was added\n")
    write(io, "# It can be reproduced/verified with md5sum if the file trimmed to the position right above [SIGNATURE]\n")
    write(io, "# The new line character preceding [SIGNATURE] BELONGS to the signature (and needs to be stripped away for recalculating/verification)\n")
    write(io, "Type $(algorithm)\n")
    write(io, "Hash $(hash_value)\n\n")
end

"""
    check_raster(seq::Sequence, raster::PulseqRaster)

"""
function check_raster(seq::Sequence, raster::PulseqRaster=DEFAULT_RASTER)
    qseq = deepcopy(seq)
    warn_count = Ref(0)
    for (bi, s) in enumerate(qseq)
        # ----- RF -----
        if is_RF_on(s)
            key = :RadiofrequencyRasterTime
            rf = qseq.RF[1, bi]
            rf.delay  = quantize_time(rf.delay,  key, getfield(raster, key), bi, "RF", "delay",  warn_count)
            rf.T      = quantize_time(rf.T,      key, getfield(raster, key), bi, "RF", "T",      warn_count; n_time_points=length(rf.A))
            rf.center = quantize_time(rf.center, key, getfield(raster, key), bi, "RF", "center", warn_count)
        end
        # ----- GR -----
        is_GR_on = [is_Gx_on(s), is_Gy_on(s), is_Gz_on(s)]
        axis = ["x", "y", "z"]
		for gi in 1:3
			if is_GR_on[gi]
                key = :GradientRasterTime
				gr = qseq.GR[gi, bi]
				gr.delay = quantize_time(gr.delay, key, getfield(raster, key), bi, "GR$(axis[gi])", "delay", warn_count)
				gr.T     = quantize_time(gr.T,     key, getfield(raster, key), bi, "GR$(axis[gi])", "T",     warn_count; n_time_points=length(gr.A))
				gr.rise  = quantize_time(gr.rise,  key, getfield(raster, key), bi, "GR$(axis[gi])", "rise",  warn_count)
				gr.fall  = quantize_time(gr.fall,  key, getfield(raster, key), bi, "GR$(axis[gi])", "fall",  warn_count)
			end
        end
        # ----- ADC -----
        adc_end = 0
        if is_ADC_on(s)
            key = :AdcRasterTime
            adc = qseq.ADC[bi]
            # Dwell time and ADC duration 
            dwell = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
            dwell = quantize_time(dwell, key, getfield(raster, key), bi, "ADC", "dwell time", warn_count)
            adc.T = adc.N == 1 ? dwell : (adc.N - 1) * dwell
            # Delay: here, there are differences between Koma and Pulseq implementations:
            # - Koma: delay = time to first sample
            # - Pulseq: delay = time to first sample - dwell/2. This is because in Pulseq the ADC event starts at a time cell edge, but samples are taken at time cell centers.
            #   Thus, we need to ensure that the Pulseq delay is a multiple of the adcRasterTime
            pulseq_delay = adc.delay - dwell/2
            if pulseq_delay < 0
                @warn "ADC delay in Koma ($(round(adc.delay*1e3, digits=4)) ms) is below dwell/2 ($(round(dwell/2*1e3, digits=4)) ms).
In Pulseq, the first ADC sample is acquired at dwell/2. Therefore, a Pulseq delay of 0 corresponds to a Koma delay of dwell/2.
This means the Koma delay must be >= dwell/2. Clamping it to this minimum value...
See https://pulseq.github.io/pulseq_shapes_and_times.pdf#page=10"
                pulseq_delay = 0
            end
            pulseq_delay = quantize_time(pulseq_delay, key, getfield(raster, key), bi, "ADC", "pulseq delay", warn_count)
            adc.delay = pulseq_delay + dwell/2
            adc_end = adc.delay + adc.T + dwell/2 # We need to add dwell/2 for the same reasons as above: the Koma ADC ends at a sample, but Pulseq ADC ends at a time cell edge.
        end
        # ----- Block -----
        key = :BlockDurationRaster
        max_event_duration = max(dur.(qseq.GR[:, bi])..., dur(qseq.RF[1, bi]), adc_end)
        block_duration = max(qseq.DUR[bi], max_event_duration)
        qseq.DUR[bi] = quantize_time(block_duration, key, getfield(raster, key), bi, "Block", "DUR", warn_count)
    end
    return qseq
end

function quantize_time(t::Number, raster_name, raster, block_id, event_key, event_element_key, warn_count; n_time_points=1)
    interval = n_time_points == 1 ? t : t / (n_time_points - 1)
    k = interval / raster
    if isapprox(k, round(k), atol=QUANT_TOL)
        return t
    else
        warn_time_quantization(warn_count, block_id, event_key, event_element_key, raster_name, raster)
        q_interval = ceil.(interval / raster) .* raster
        return n_time_points == 1 ? q_interval : q_interval * (n_time_points - 1)
    end
end

function quantize_time(t::Array, raster_name, raster, block_id, event_key, event_element_key, warn_count; n_time_points=1)
    k = t / raster
    if all(isapprox(k_i, round(k_i), atol=QUANT_TOL) for k_i in k)
        return t
    else
        warn_time_quantization(warn_count, block_id, event_key, event_element_key, raster_name, raster)
        return ceil.(t / raster) .* raster
    end
end

function warn_time_quantization(warn_count, block_id, event_key, event_element_key, raster_name, raster)
    if warn_count[] < 10
        @warn "Block $block_id: Event $event_key: Element $event_element_key:
Time is not a multiple of $(string(raster_name)) ($(raster * 1e6) μs). Quantizing it..."
    elseif warn_count[] == 10
        @warn "Additional time quantization warnings occurred; detailed logs are capped at 10."
    end
    warn_count[] += 1
end

"""
    write_seq(seq, filename)

Writes a Sequence struct to a Pulseq file with `.seq` extension.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `filename`: (`::String`) absolute or relative path of the sequence file `.seq`

# Examples
```julia-repl
julia> seq = Sequence()
julia> seq += RF(2e-6, 3e-3)
julia> seq += ADC(100, 100e-3)
julia> write_seq(seq, "fid.seq")
```
"""
function write_seq(
    seq::Sequence, filename::String;
    sys::Scanner=Scanner(),
    signatureAlgorithm::AbstractString = "md5"
)
    # 1. Check scanner constraints. If the sequence is not compliant, the function will throw an error.
    @info "Checking scanner constraints..." B0_max = sys.B0 B1_max = sys.B1 G_max = sys.Gmax S_max = sys.Smax ADC_Δt = sys.ADC_Δt
    check_scanner_constraints(seq, sys)
    # 2. Create the Pulseq Raster 
    raster = PulseqRaster(seq, sys)
    # 3. Check raster. If the sequence is not compliant, the function will throw a warning and return the sequence with the correct raster.
    @info "Checking Pulseq raster..." BlockDurationRaster = raster.BlockDurationRaster GradientRasterTime = raster.GradientRasterTime RadiofrequencyRasterTime = raster.RadiofrequencyRasterTime AdcRasterTime = raster.AdcRasterTime
    seq = check_raster(seq, raster)
    @info "Saving sequence to $(basename(filename)) ..."
    # 4. Collect the pulseq assets.
    blocks, event_libraries = collect_pulseq_assets(seq, raster)
    buffer = IOBuffer()
    emit_pulseq(buffer, blocks, event_libraries)
    payload = take!(buffer)
    signature_hash = supported_signature_digest(signatureAlgorithm, payload)
    open(filename, "w") do io
        write(io, payload)
        write(io, '\n')
        emit_signature_section!(io, signatureAlgorithm, signature_hash)
    end
    return nothing
end
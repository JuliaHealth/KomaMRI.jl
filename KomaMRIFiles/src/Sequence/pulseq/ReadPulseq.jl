parse_keyed_value(line, key, T) = begin
    fields = split(line)
    length(fields) == 2 && fields[1] == key || error("Expected `$key <value>`, got `$line`.")
    return parse(T, fields[2])
end

try_parse_keyed_value(line, key, T) = begin
    fields = split(line)
    length(fields) == 2 && fields[1] == key || return nothing
    return parse(T, fields[2])
end

"""
read_version Read the [VERSION] section of a sequence file.
   defs=read_version(fid) Read Pulseq version from file
   identifier of an open MR sequence file and return it
"""
function read_version(io; verbose=true)
    pulseq_version = VersionNumber(
        parse_keyed_value(readline(io), "major", Int),
        parse_keyed_value(readline(io), "minor", Int),
        parse_keyed_value(readline(io), "revision", Int),
    )
    @assert pulseq_version.major == 1 "Unsupported version_major $(pulseq_version.major)"
    if     pulseq_version < v"1.2.0"
        @error "Unsupported Pulseq $(pulseq_version), only file format revision 1.2.0 and above are supported"
    elseif pulseq_version < v"1.3.1"
        @warn "Loading older Pulseq $(pulseq_version); some code may not function as expected"
    elseif pulseq_version >= v"1.5.0" && verbose
        @info "Pulseq $(pulseq_version) is supported, but Soft Delay and RF Shimming extensions are not yet included\n(see https://github.com/JuliaHealth/KomaMRI.jl/issues/714)" maxlog=1
    end
    pulseq_version
end

"""
read_definitions Read the [DEFINITIONS] section of a sequence file.
   defs=read_definitions(fid) Read user definitions from file
   identifier of an open MR sequence file and return a map of
   key/value entries.
"""
function read_definitions(io)
    entries = Pair{String,String}[]
    block_duration_raster = DEFAULT_RASTER.BlockDurationRaster
    gradient_raster_time = DEFAULT_RASTER.GradientRasterTime
    radiofrequency_raster_time = DEFAULT_RASTER.RadiofrequencyRasterTime
    adc_raster_time = DEFAULT_RASTER.AdcRasterTime
    required_extensions = String[]
    while true
        line = readline(io)
        fields = split(line)
        isempty(fields) && break
        key = String(fields[1])
        value = join(fields[2:end], " ")
        push!(entries, key => value)
        if key == "BlockDurationRaster"
            block_duration_raster = parse(Float64, fields[2])
        elseif key == "GradientRasterTime"
            gradient_raster_time = parse(Float64, fields[2])
        elseif key == "RadiofrequencyRasterTime"
            radiofrequency_raster_time = parse(Float64, fields[2])
        elseif key == "AdcRasterTime"
            adc_raster_time = parse(Float64, fields[2])
        elseif key == "RequiredExtensions"
            required_extensions = String.(fields[2:end])
        end
    end
    return PulseqDefinitions(
        entries,
        block_duration_raster,
        gradient_raster_time,
        radiofrequency_raster_time,
        adc_raster_time,
        required_extensions,
    )
end

function parsed_definition_value(key, value)
    fields = split(value)
    parsed = map(fields) do field
        something(tryparse(Float64, field), field)
    end
    return length(parsed) == 1 && key != "RequiredExtensions" ? parsed[1] : parsed
end

function definitions_dict(definitions)
    dict = Dict{String,Any}()
    for (key, value) in definitions.entries
        dict[key] = parsed_definition_value(key, value)
    end
    dict["BlockDurationRaster"] = definitions.block_duration_raster
    dict["GradientRasterTime"] = definitions.gradient_raster_time
    dict["RadiofrequencyRasterTime"] = definitions.radiofrequency_raster_time
    dict["AdcRasterTime"] = definitions.adc_raster_time
    isempty(definitions.required_extensions) || (dict["RequiredExtensions"] = copy(definitions.required_extensions))
    return dict
end

"""
read_signature Read the [SIGNATURE] section of a sequence file.
   signature=read_signature(fid) Read user signature from file
   identifier of an open MR sequence file and return a NamedTuple with signature metadata.
"""
function read_signature(io)
    signature_type = nothing
    signature_hash = nothing
    while true
        line = readline(io)
        fields = split(line)
        isempty(fields) && break
        firstchar = first(fields[1])
        firstchar == '#' && continue
        key = fields[1]
        value = join(fields[2:end], " ")
        if key == "Type"
            signature_type = strip(value)
        elseif key == "Hash"
            signature_hash = lowercase(replace(strip(value), " " => ""))
        end
    end
    has_hash = !isnothing(signature_type) && !isnothing(signature_hash)
    return has_hash ? (type = signature_type, hash = signature_hash) : nothing
end


"""
read_blocks Read the [BLOCKS] section of a sequence file.
   library=read_blocks(fid, blockDurationRaster, pulseq_version)
   Read blocks from file identifier of an open MR sequence file and
   returns block IDs, durations and delay IDs as arrays or vectors.

"""
function read_blocks(io, block_duration_raster, pulseq_version)
    n_block_events = pulseq_version <= v"1.2.1" ? 7 : 8
    blocks = Int[]
    blockDurations = Float64[]
    delayIDs_tmp = Int[]
    block_events = Vector{Int}(undef, n_block_events)
    num_blocks = 0
    while true
        line = readline(io)
        isempty(line) && break
        n_parsed = parse_block_events!(block_events, line)
        if n_parsed != n_block_events
            error("Expected $n_block_events events but got $n_parsed.")
        end
        if block_events[1] != 0
            if pulseq_version <= v"1.2.1"
                # discard id and duration (duration treated below)
                append!(blocks, @view(block_events[3:end]))
                push!(blocks, 0) # append 0 for version compatibility (see assert below)
            else
                append!(blocks, @view(block_events[3:end]))
            end
            if pulseq_version >= v"1.4.0" # Explicit block duration (in units of blockDurationRaster)
                duration = block_events[2] * block_duration_raster
                push!(blockDurations, duration)
            else # Implicit block duration from delay ID
                push!(delayIDs_tmp, block_events[2])
            end
            num_blocks += 1
        end
    end
    reshaped_blocks = reshape(blocks, :, num_blocks)
    # we need 6 vals because we drop IDs and durations
    @assert size(reshaped_blocks, 1) == 6 "unexpected number of fields per block"
    return reshaped_blocks, blockDurations, delayIDs_tmp
end

function parse_block_events!(dest, line)
    n = 0
    value = 0
    sign = 1
    in_number = false
    for byte in codeunits(line)
        if byte == UInt8(' ') || byte == UInt8('\t') || byte == UInt8('\r')
            if in_number
                n += 1
                dest[n] = sign * value
                value = 0
                sign = 1
                in_number = false
            end
        elseif byte == UInt8('-')
            sign = -1
            in_number = true
        else
            value = 10 * value + Int(byte - UInt8('0'))
            in_number = true
        end
    end
    if in_number
        n += 1
        dest[n] = sign * value
    end
    return n
end

apply_scale(scale::Number, value::Number) = scale * value
apply_scale(scale, value) = value
parse_event_field(token, field) =
    token == "%i" ? parse(Int, field) :
    token == "%f" ? parse(Float64, field) :
    token == "%c" ? only(field) :
    String(field)

"""
    events = read_events(io, scale; field_format="%i "*"%f "^(length(scale)), event_library)

Read an event section of a sequence file.

# Arguments
- `io`: (`::IO`) input file stream
- `scale`: tuple of scale factors applied to parsed fields

# Keywords
- `field_format`: (`::String`) Pulseq field format tokens (`%i`, `%f`, `%c`, `%s`)
- `event_library`: output dictionary for the parsed event library
- `constructor`: constructor used to normalize each parsed event row

# Returns
- `events`: parsed event library
"""
function read_events(io, scales; field_format="%i " * "%f " ^ length(scales), event_library, constructor=identity)
    field_formats = split(field_format)
    event_length = length(scales) + 1
    while true
        line = readline(io)
        isempty(line) && break
        fields = split(line)
        length(fields) == event_length || break #Break if not all values read
        id = parse(Int, fields[1])
        parsed = ntuple(length(scales)) do i
            apply_scale(scales[i], parse_event_field(field_formats[i + 1], fields[i + 1]))
        end
        event_library[id] = constructor(parsed)
    end
    return event_library
end

"""
    read_extensions(io, ext_string, ext_type, ext_id, extensionTypeLibrary, extensionSpecLibrary, required_extensions)

This function will read the extension specifications (for the example below, the lines 
after "extension TRRIGERS 1") and add them to the `extensionSpecLibrary` dictionary.
It will also add the extension type to the `extensionTypeLibrary` dictionary.

# Example
Pulseq example file:
```
[EXTENSIONS]
1 1 1 0
2 1 2 0

extension TRIGGERS 1
1 2 1 0 2000
2 1 3 500 100
```
"""
function read_extensions(io, ext_string, ext_type::Type{<:Extension}, ext_id, extensionTypeLibrary, extensionSpecLibrary, required_extensions)
    extensionTypeLibrary[ext_id] = ext_type
    raw_specs = read_events(
        io,
        Tuple(get_scale(ext_type));
        field_format="%i " * get_pulseq_format(ext_type),
        event_library=Dict{Int, Tuple}(),
    )
    extensionSpecLibrary[ext_id] = Dict(id => ext_type(values...) for (id, values) in raw_specs)
end
function read_extensions(io, ext_string, ext_type, ext_id, extensionTypeLibrary, extensionSpecLibrary, required_extensions)
    if ext_string in required_extensions
        @warn "Ignoring unsupported required extension: $ext_string" RequiredExtensions = required_extensions
    else
        @warn "Ignoring unsupported extension: $ext_string"
    end
    while true # Skip the extension specifications
        line = readline(io)
        isempty(line) && break
    end
end

"""
read_shapes Read the [SHAPES] section of a sequence file.
   library=read_shapes(fid) Read shapes from file identifier of an
   open MR sequence file and return a library of shapes.
"""
function read_shapes(io, forceConvertUncompressed)
    shapeLibrary = ShapeLibrary()
    while true #Reading shapes
        eof(io) && break
        mark(io) # Mark the position before reading the line
        line = readline(io)
        if isempty(line) unmark(io); continue end
        id = try_parse_keyed_value(line, "shape_id", Int)
        if isnothing(id)
            reset(io)
            break
        end
        unmark(io) # Unmark the position after reading the line
        num_samples = parse_keyed_value(readline(io), "num_samples", Int)
        shape = Float64[]
        sizehint!(shape, num_samples)
        while true #Reading shape data
            data_point = tryparse(Float64, readline(io))
            !isnothing(data_point) || break #Break if no sample
            push!(shape, data_point)
        end
        # check if conversion is needed: in v1.4.x we use length(data)==num_samples
        # as a marker for the uncompressed (stored) data. In older versions this condition could occur by chance
        if forceConvertUncompressed && length(shape)==num_samples
            num_samples, shape = compress_shape(decompress_shape(num_samples, shape; forceDecompression=true))
        end
        data = (num_samples, shape)
        shapeLibrary[id] = data
    end
    return shapeLibrary
end

"""
    fix_first_last_grads!(blockEvents, blockDurations, eventLibraries)

Updates the `eventLibraries` dictionary with new first and last points for gradients.

# Notes:
- This function is "replicating" the following MATLAB code:
https://github.com/pulseq/pulseq/blob/v1.5.1/matlab/%2Bmr/%40Sequence/read.m#L325-L413
- We are updating the `gradLibrary` entries with the new first and last points, making them compatible with the v1.5.x format.
"""
function fix_first_last_grads!(blockEvents, blockDurations, eventLibraries)
    # Add first and last Pulseq points
    grad_prev_last = zeros(3)
    grad_durations = Dict{Int,Float64}()
    for (iB, eventIDs) in enumerate(eachcol(blockEvents))
        block_duration = blockDurations[iB]
        for iG in eachindex(grad_prev_last)
            g_id = eventIDs[1 + iG]
            g_id > 0 || continue
            update_first_last!(grad_prev_last, iG, eventLibraries.grad_library[g_id], block_duration, eventLibraries, grad_durations, g_id)
        end
    end
end

function fix_legacy_trapezoids!(grad_library, grad_raster_time)
    for (i, grad) in grad_library
        grad isa PulseqTrapGradEvent || continue
        fix_rise = grad.rise == 0 && abs(grad.amplitude) == 0 && grad.flat > 0
        flat = fix_rise ? grad.flat - grad_raster_time : grad.flat
        fix_delay = grad.delay == 0 && abs(grad.amplitude) == 0 && flat > 0
        (fix_rise || fix_delay) || continue
        grad_library[i] = PulseqTrapGradEvent(
            grad.amplitude,
            fix_rise ? grad_raster_time : grad.rise,
            fix_delay ? flat - grad_raster_time : flat,
            grad.fall,
            fix_delay ? grad_raster_time : grad.delay,
        )
    end
end

# Decoded runtime events used to materialize concrete Sequence storage.
struct PulseqDecodedLibraries{G,R,H,A,E}
    gradients::G
    rfs::R
    rf_half_steps::H
    adcs::A
    extensions::E
end

const PulseqDecodedExtensions = Vector{Extension}
const PULSEQ_EMPTY_GRAD = Grad(0.0, 0.0)
const PULSEQ_EMPTY_RF = RF(0.0, 0.0)
const PULSEQ_EMPTY_ADC = ADC(0, 0.0)
const PULSEQ_EMPTY_EXTENSIONS = PulseqDecodedExtensions()

@inline decoded_or_empty(entries, id::Int, empty) = iszero(id) ? copy(empty) : copy(entries[id])
decoded_event_library(f, n, empty) = iszero(n) ? typeof(empty)[] : [f(id) for id in 1:n]
function decoded_small_union_eltype(entries, empty)
    types = DataType[typeof(empty)]
    for entry in entries
        T = typeof(entry)
        T in types || push!(types, T)
    end
    return length(types) == 1 ? only(types) : Core.apply_type(Union, types...)
end

function init_legacy_block_durations!(blockDurations, blockEvents, delayIDs_tmp, eventLibraries, decodedLibraries)
    resize!(blockDurations, size(blockEvents, 2))
    for i in axes(blockEvents, 2)
        delayID = delayIDs_tmp[i]
        delay = delayID > 0 ? eventLibraries.tmp_delay_library[delayID] : 0.0
        Gx, Gy, Gz, rf, add_half_Δt_rf, adc, _ = decoded_block(
            decodedLibraries,
            blockEvents[1, i],
            blockEvents[2, i],
            blockEvents[3, i],
            blockEvents[4, i],
            blockEvents[5, i],
            blockEvents[6, i],
        )
        blockDurations[i] = max(
            delay,
            dur(Gx),
            dur(Gy),
            dur(Gz),
            dur(rf) + add_half_Δt_rf * eventLibraries.definitions.radiofrequency_raster_time / 2,
            dur(adc),
        )
    end
end

max_pulseq_id(dict) = isempty(dict) ? 0 : maximum(keys(dict))

function decode_pulseq_libraries(eventLibraries)
    n_grads = max_pulseq_id(eventLibraries.grad_library)
    grad_library = decoded_event_library(n_grads, PULSEQ_EMPTY_GRAD) do id
        get_Grad(eventLibraries.grad_library, eventLibraries.shape_library, eventLibraries.definitions.gradient_raster_time, id)
    end

    n_rfs = max_pulseq_id(eventLibraries.rf_library)
    rf_library = decoded_event_library(n_rfs, PULSEQ_EMPTY_RF) do id
        get_RF(eventLibraries.rf_library, eventLibraries.shape_library, eventLibraries.definitions.radiofrequency_raster_time, id)
    end
    rf_half_steps = [pulseq_rf_adds_half_step(eventLibraries.rf_library[id]) for id in 1:n_rfs]

    n_adcs = max_pulseq_id(eventLibraries.adc_library)
    adc_library = Vector{ADC}(undef, n_adcs)
    for id in 1:n_adcs
        adc_library[id] = get_ADC(eventLibraries.adc_library, id)
    end

    n_extensions = max_pulseq_id(eventLibraries.extension_instance_library)
    extension_library = Vector{PulseqDecodedExtensions}(undef, n_extensions)
    for id in 1:n_extensions
        extension_library[id] = get_extension(
            eventLibraries.extension_instance_library,
            eventLibraries.extension_type_library,
            eventLibraries.extension_spec_library,
            id,
        )
    end

    return PulseqDecodedLibraries(
        grad_library,
        rf_library,
        rf_half_steps,
        adc_library,
        extension_library,
    )
end

function update_first_last!(grad_prev_last, iG, ::PulseqTrapGradEvent, block_duration, event_libraries, grad_durations, g_id)
    grad_prev_last[iG] = 0.0
    return nothing
end

function update_first_last!(grad_prev_last, iG, g::PulseqArbGradEvent, block_duration, event_libraries, grad_durations, g_id)
    grad_duration = get!(grad_durations, g_id) do
        updated = update_first_last_grad(g, grad_prev_last[iG], event_libraries.shape_library)
        event_libraries.grad_library[g_id] = updated
        pulseq_event_duration(updated, event_libraries.shape_library, event_libraries.definitions.gradient_raster_time)
    end
    g = event_libraries.grad_library[g_id]
    grad_prev_last[iG] = grad_duration + eps(Float64) < block_duration ? 0 : g.last
    return nothing
end

update_first_last_grad(g::PulseqArbGradEvent, first, shape_library) = PulseqArbGradEvent(
    g.amplitude,
    first,
    pulseq_last_gradient_sample(first, g, shape_library),
    g.amp_shape_id,
    g.time_shape_id,
    g.delay,
)

function pulseq_last_gradient_sample(first, g::PulseqArbGradEvent, shape_library)
    waveform = g.amplitude * decompress_shape(shape_library[g.amp_shape_id]...)
    if g.time_shape_id != 0
        return waveform[end]
    end
    last = first
    for sample in waveform
        last = 2 * sample - last
    end
    return last
end

pulseq_event_duration(g::PulseqTrapGradEvent, shape_library, Δt_gr) = g.delay + g.rise + g.flat + g.fall

function pulseq_event_duration(g::PulseqArbGradEvent, shape_library, Δt_gr)
    num_samples = shape_library[g.amp_shape_id][1]
    if g.time_shape_id == 0
        return g.delay + num_samples * Δt_gr
    elseif g.time_shape_id == -1
        return g.delay + num_samples * (Δt_gr / 2)
    else
        tt = decompress_shape(shape_library[g.time_shape_id]...)
        return g.delay + tt[end] * Δt_gr
    end
end

@inline function decoded_block(decodedLibraries::PulseqDecodedLibraries, irf::Int, igx::Int, igy::Int, igz::Int, iadc::Int, iext::Int)
    Gx = decoded_or_empty(decodedLibraries.gradients, igx, PULSEQ_EMPTY_GRAD)
    Gy = decoded_or_empty(decodedLibraries.gradients, igy, PULSEQ_EMPTY_GRAD)
    Gz = decoded_or_empty(decodedLibraries.gradients, igz, PULSEQ_EMPTY_GRAD)
    rf = decoded_or_empty(decodedLibraries.rfs, irf, PULSEQ_EMPTY_RF)
    adc = decoded_or_empty(decodedLibraries.adcs, iadc, PULSEQ_EMPTY_ADC)
    ext = decoded_or_empty(decodedLibraries.extensions, iext, PULSEQ_EMPTY_EXTENSIONS)
    add_half_Δt_rf = irf != 0 && decodedLibraries.rf_half_steps[irf]
    return Gx, Gy, Gz, rf, add_half_Δt_rf, adc, ext
end

# Decode dense block-id tables into concrete runtime events and build the Sequence.
function get_seq_from_blocks(blocks::Vector{PulseqBlockEventIDs}, definitions, decodedLibraries)
    num_blocks = length(blocks)
    # Preserve Pulseq event families without falling back to abstract Matrix{Grad}/Matrix{RF}.
    GR = Matrix{decoded_small_union_eltype(decodedLibraries.gradients, PULSEQ_EMPTY_GRAD)}(undef, 3, num_blocks)
    RFs = Matrix{decoded_small_union_eltype(decodedLibraries.rfs, PULSEQ_EMPTY_RF)}(undef, 1, num_blocks)
    ADCs = Vector{ADC}(undef, num_blocks)
    DUR = Vector{Float64}(undef, num_blocks)
    EXTs = Vector{PulseqDecodedExtensions}(undef, num_blocks)

    if Threads.nthreads() > 1
        Threads.@threads for i in eachindex(blocks)
            fill_decoded_block!(GR, RFs, ADCs, DUR, EXTs, blocks, definitions, decodedLibraries, i)
        end
    else
        for i in eachindex(blocks)
            fill_decoded_block!(GR, RFs, ADCs, DUR, EXTs, blocks, definitions, decodedLibraries, i)
        end
    end

    return KomaMRIBase.Sequence(GR, RFs, ADCs, DUR, EXTs, Dict{String,Any}())
end

@inline function fill_decoded_block!(GR, RFs, ADCs, DUR, EXTs, blocks, definitions, decodedLibraries, i)
    block = blocks[i]
    Gx, Gy, Gz, rf, _, adc, ext = decoded_block(
        decodedLibraries,
        block.rf_id,
        block.gx_id,
        block.gy_id,
        block.gz_id,
        block.adc_id,
        block.ext_id,
    )
    GR[1, i] = Gx
    GR[2, i] = Gy
    GR[3, i] = Gz
    RFs[1, i] = rf
    ADCs[i] = adc
    DUR[i] = block.duration_ticks * definitions.block_duration_raster
    EXTs[i] = ext
    return nothing
end

"""
    PulseqSequenceData

Pulseq-native sequence representation used by [`read_seq_data`](@ref) and
[`write_seq_data`](@ref). It stores block event ids, event libraries, definitions,
the Pulseq version, and the optional signature without materializing repeated
Koma `RF`, `Grad`, and `ADC` objects.

Event amplitudes are stored in SI units (RF: T, gradients: T/m). Serialized time
fields are normalized from μs/ns to seconds.
"""
struct PulseqSequenceData
    blocks::Vector{PulseqBlockEventIDs}
    libraries::PulseqFileEventLibraries
    pulseq_version::VersionNumber
    signature::Union{Nothing,PulseqSignature}
end

# Section-by-section parse accumulator before version-specific normalization.
mutable struct PulseqParsedFile
    pulseq_version::VersionNumber
    definitions::PulseqDefinitions
    signature::Union{Nothing,PulseqSignature}
    blockEvents::Matrix{Int}
    blockDurations::Vector{Float64}
    delayIDs_tmp::Vector{Int}
    gradLibrary::Dict{Int,PulseqGradEvent}
    rfLibrary::Dict{Int,PulseqRFEvent}
    adcLibrary::Dict{Int,PulseqADCEvent}
    tmp_delayLibrary::Dict{Int,Float64}
    shapeLibrary::ShapeLibrary
    extensionInstanceLibrary::Dict{Int,PulseqExtensionInstanceEvent}
    extensionTypeLibrary::Dict{Int,Type{<:Extension}}
    extensionSpecLibrary::Dict{Int,Dict{Int,Extension}}
end

PulseqParsedFile() = PulseqParsedFile(
    v"0.0.0",
    PulseqDefinitions(),
    nothing,
    Matrix{Int}(undef, 0, 0),
    Float64[],
    Int[],
    Dict{Int,PulseqGradEvent}(),
    Dict{Int,PulseqRFEvent}(),
    Dict{Int,PulseqADCEvent}(),
    Dict{Int,Float64}(),
    ShapeLibrary(),
    Dict{Int,PulseqExtensionInstanceEvent}(),
    Dict{Int,Type{<:Extension}}(),
    Dict{Int,Dict{Int,Extension}}(),
)

"""
    parsed = parse_pulseq_file(io)

Parse one Pulseq file stream into the intermediate `PulseqParsedFile` representation.
Used by [`read_seq_data`](@ref).
"""
function parse_pulseq_file(io; verbose=true)
    parsed = PulseqParsedFile()
    while !eof(io)
        section = readline(io)
        if isempty(section) || section[1] == '#'
            continue
        elseif section == "[DEFINITIONS]"
            parsed.definitions = read_definitions(io)
        elseif section == "[VERSION]"
            parsed.pulseq_version = read_version(io; verbose)
        elseif section == "[BLOCKS]"
            if parsed.pulseq_version == v"0.0.0"
                @error "Pulseq file MUST include [VERSION] section prior to [BLOCKS] section"
            end
            parsed.blockEvents, parsed.blockDurations, parsed.delayIDs_tmp =
                read_blocks(io, parsed.definitions.block_duration_raster, parsed.pulseq_version)
        elseif section == "[RF]"
            if parsed.pulseq_version >= v"1.5.0"
                parsed.rfLibrary = read_events(io, (1 / γ, 1, 1, 1, 1e-6, 1e-6, 1, 1, 1, 1, 1); field_format="%i " * "%f " ^ 10 * "%c ", event_library=parsed.rfLibrary, constructor=PulseqRFEvent)
            elseif parsed.pulseq_version >= v"1.4.0"
                parsed.rfLibrary = read_events(io, (1 / γ, 1, 1, 1, 1e-6, 1, 1); event_library=parsed.rfLibrary, constructor=PulseqRFEvent)
            else
                parsed.rfLibrary = read_events(io, (1 / γ, 1, 1, 1e-6, 1, 1); event_library=parsed.rfLibrary, constructor=PulseqRFEvent)
            end
        elseif section == "[GRADIENTS]"
            if parsed.pulseq_version >= v"1.5.0"
                parsed.gradLibrary = read_events(io, (1 / γ, 1 / γ, 1 / γ, 1, 1, 1e-6); event_library=parsed.gradLibrary, constructor=PulseqArbGradEvent)
            elseif parsed.pulseq_version >= v"1.4.0"
                parsed.gradLibrary = read_events(io, (1 / γ, 1, 1, 1e-6); event_library=parsed.gradLibrary, constructor=PulseqArbGradEvent)
            else
                parsed.gradLibrary = read_events(io, (1 / γ, 1, 1e-6); event_library=parsed.gradLibrary, constructor=PulseqArbGradEvent)
            end
        elseif section == "[TRAP]"
            parsed.gradLibrary = read_events(io, (1 / γ, 1e-6, 1e-6, 1e-6, 1e-6); event_library=parsed.gradLibrary, constructor=PulseqTrapGradEvent)
        elseif section == "[ADC]"
            if parsed.pulseq_version >= v"1.5.0"
                parsed.adcLibrary = read_events(io, (1, 1e-9, 1e-6, 1, 1, 1, 1, 1); event_library=parsed.adcLibrary, constructor=PulseqADCEvent)
            else
                parsed.adcLibrary = read_events(io, (1, 1e-9, 1e-6, 1, 1); event_library=parsed.adcLibrary, constructor=PulseqADCEvent)
            end
        elseif section == "[DELAYS]"
            if parsed.pulseq_version >= v"1.4.0"
                @error "Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section"
            end
            parsed.tmp_delayLibrary = read_events(io, (1e-6,); event_library=parsed.tmp_delayLibrary, constructor=first)
        elseif section == "[SHAPES]"
            parsed.shapeLibrary = read_shapes(io, parsed.pulseq_version.major == 1 && parsed.pulseq_version.minor < 4)
        elseif section == "[EXTENSIONS]"
            parsed.extensionInstanceLibrary = read_events(io, (1, 1, 1); field_format="%i " * "%i " ^ 3, event_library=parsed.extensionInstanceLibrary, constructor=PulseqExtensionInstanceEvent)
        elseif section == "[SIGNATURE]"
            parsed.signature = read_signature(io)
        elseif startswith(section, "extension")
            ext = section[11:end]
            fields = split(ext)
            ext_string = String(fields[1])
            ext_type = get_EXT_type_from_symbol(Val(Symbol(ext_string)))
            ext_id = parse(Int, fields[2])
            read_extensions(io, ext_string, ext_type, ext_id, parsed.extensionTypeLibrary, parsed.extensionSpecLibrary, parsed.definitions.required_extensions)
        else
            @error "Unknown section code: $section"
        end
    end
    return parsed
end

function pulseq_event_libraries(parsed::PulseqParsedFile)
    return PulseqFileEventLibraries(
        parsed.gradLibrary,
        parsed.rfLibrary,
        parsed.adcLibrary,
        parsed.tmp_delayLibrary,
        parsed.shapeLibrary,
        parsed.extensionInstanceLibrary,
        parsed.extensionTypeLibrary,
        parsed.extensionSpecLibrary,
        parsed.definitions,
    )
end

function pulseq_blocks(blockEvents, blockDurations, definitions)
    blocks = Vector{PulseqBlockEventIDs}(undef, size(blockEvents, 2))
    for i in axes(blockEvents, 2)
        blocks[i] = PulseqBlockEventIDs(
            i,
            pulseq_block_duration_ticks(blockDurations[i], definitions.block_duration_raster),
            blockEvents[1, i],
            blockEvents[2, i],
            blockEvents[3, i],
            blockEvents[4, i],
            blockEvents[5, i],
            blockEvents[6, i],
        )
    end
    return blocks
end

function pulseq_sequence_data(parsed::PulseqParsedFile; filename=nothing, verify_signature=false)
    pulseq_version = parsed.pulseq_version
    signature = parsed.signature
    blockEvents = parsed.blockEvents
    blockDurations = parsed.blockDurations
    delayIDs_tmp = parsed.delayIDs_tmp
    gradLibrary = parsed.gradLibrary

    # fix trapezoidal gradients imported from older versions
    if pulseq_version < v"1.4.0"
        fix_legacy_trapezoids!(gradLibrary, parsed.definitions.gradient_raster_time)
    end
    verify_signature && !isnothing(filename) && verify_signature!(filename, signature; pulseq_version=pulseq_version)
    eventLibraries = pulseq_event_libraries(parsed)

    if pulseq_version < v"1.4.0"
        decodedLibraries = decode_pulseq_libraries(eventLibraries)
        init_legacy_block_durations!(blockDurations, blockEvents, delayIDs_tmp, eventLibraries, decodedLibraries)
    end
    # Add first and last points for gradients #320 for version <= 1.4.2
    if pulseq_version < v"1.5.0"
        fix_first_last_grads!(blockEvents, blockDurations, eventLibraries)
    end

    blocks = pulseq_blocks(blockEvents, blockDurations, parsed.definitions)
    return PulseqSequenceData(blocks, eventLibraries, pulseq_version, signature)
end

"""
    data = read_seq_data(filename)

Read a Pulseq file into Pulseq's native block/event-library representation without
materializing repeated Koma `Sequence` events.
"""
function read_seq_data(filename::AbstractString; verify_signature=false, verbose=true)
    parsed = open(filename) do io
        parse_pulseq_file(io; verbose)
    end
    return pulseq_sequence_data(parsed; filename, verify_signature)
end

function read_seq_data(io::IO; verify_signature=false, verbose=true)
    parsed = parse_pulseq_file(io; verbose)
    return pulseq_sequence_data(parsed; verify_signature)
end

function sequence_from_pulseq_data(data::PulseqSequenceData; filename=nothing)
    decodedLibraries = decode_pulseq_libraries(data.libraries)
    seq = get_seq_from_blocks(data.blocks, data.libraries.definitions, decodedLibraries)

    # Final details
    seq.DEF = KomaMRIBase._sequence_def_from_pulseq(definitions_dict(data.libraries.definitions))
    # Koma specific details for reconstrucion
    if !isnothing(filename)
        seq.DEF["FileName"] = basename(filename)
        seq.DEF["PulseqVersion"] = data.pulseq_version
        seq.DEF["signature"] = data.signature
    end
    # Guessing recon dimensions
    seq.DEF["Nx"] = trunc(Int64, get(seq.DEF, "Nx", maximum(adc.N for adc = seq.ADC)))
    seq.DEF["Nz"] = trunc(Int64, get(seq.DEF, "Nz", length(unique(seq.RF.Δf))))
    seq.DEF["Ny"] = trunc(Int64, get(seq.DEF, "Ny", count(is_ADC_on, seq.ADC) ÷ seq.DEF["Nz"]))
    #Koma sequence
    return seq
end

KomaMRIBase.Sequence(data::PulseqSequenceData; filename=nothing) =
    sequence_from_pulseq_data(data; filename)

"""
    seq = read_seq(filename)

Returns the Sequence struct from a Pulseq file with `.seq` extension.

# Arguments
- `filename`: (`::String`) absolute or relative path of the sequence file `.seq`

# Keywords
- `verify_signature`: (`::Bool`, `=false`) verify the optional Pulseq `[SIGNATURE]` hash
  while loading
- `apply_rotations`: (`::Bool`, `=true`) apply Pulseq `ROTATIONS` extensions to
  the gradients while keeping the extension events in `seq.EXT`
- `verbose`: (`::Bool`, `=true`) show informational loading messages

# Returns
- `seq`: (`::Sequence`) Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_seq(seq)
```
"""
function read_seq(filename; verify_signature=false, apply_rotations=true, verbose=true)
    verbose && @info "Loading sequence $(basename(filename)) ..."
    data = read_seq_data(filename; verify_signature, verbose)
    seq = KomaMRIBase.Sequence(data; filename)
    has_rotations = any(==(KomaMRIBase.QuaternionRot), values(data.libraries.extension_type_library))
    return apply_rotations && has_rotations ?
        KomaMRIBase._apply_rotations_to_owned_sequence(seq) :
        seq
end

#To Sequence
"""
    grad = get_Grad(gradLibrary, shapeLibrary, Δt_gr, i)

Decode one Pulseq gradient event into a runtime `Grad`. Used by
[`decode_pulseq_libraries`](@ref).

# Arguments
- `gradLibrary`: (`::Dict{Int, PulseqGradEvent}`) the "gradLibrary" dictionary
- `shapeLibrary`: (`::ShapeLibrary`) the "shapeLibrary" dictionary
- `Δt_gr`: (`::Float64`, `[s]`) gradient raster time
- `i`: (`::Int64`) index of the axis in the block event

# Returns
- `grad`: (::Grad) Gradient struct
"""
function get_Grad(gradLibrary, shapeLibrary, Δt_gr, i)
    return get_Grad(gradLibrary[i], shapeLibrary, Δt_gr)
end

get_Grad(grad::PulseqTrapGradEvent, shapeLibrary, Δt_gr) =
    Grad(grad.amplitude, grad.flat, grad.rise, grad.fall, grad.delay, 0.0, 0.0)

function get_Grad(grad::PulseqArbGradEvent, shapeLibrary, Δt_gr)
    amplitude = grad.amplitude
    first_grads = grad.first
    last_grads = grad.last
    amp_shape_id = grad.amp_shape_id
    time_shape_id = grad.time_shape_id
    delay = grad.delay
    #Amplitude
    gA = amplitude * decompress_shape(shapeLibrary[amp_shape_id]...)
    n_gr = length(gA) - 1
    #Creating timings
    if time_shape_id == 0 #no time waveform. Default time raster
        gT = n_gr * Δt_gr
        rise, fall = Δt_gr/2, Δt_gr/2
    elseif time_shape_id == -1 #New in pulseq 1.5.x: no time waveform. 1/2 of the default time raster
        gT = n_gr * Δt_gr / 2
        rise, fall = Δt_gr/2, Δt_gr/2
    else
        gt = decompress_shape(shapeLibrary[time_shape_id]...)
        gt[1] == 0 || @warn "Gradient time shape $time_shape_id starting at a non-zero value $(gt[1]). This is not recommended and may not be supported properly\n (see https://github.com/pulseq/pulseq/issues/188#issuecomment-3541588756) " maxlog=1
        delay += gt[1] * Δt_gr # offset due to the shape starting at a non-zero value. This case 
        gT = diff(gt) * Δt_gr
        rise, fall = 0.0, 0.0
    end
    return Grad(gA, gT, rise, fall, delay, first_grads, last_grads)
end

"""
    rf = get_RF(rfLibrary, shapeLibrary, Δt_rf, i)

Decode one Pulseq RF event into a runtime `RF`. Used by
[`decode_pulseq_libraries`](@ref).

# Arguments
- `rfLibrary`: (`::Dict{Int, PulseqRFEvent}`) the "rfLibrary" dictionary
- `shapeLibrary`: (`::ShapeLibrary`) the "shapeLibrary" dictionary
- `Δt_rf`: (`::Float64`, `[s]`) RF raster time
- `i`: (`::Int64`) index of the RF in the block event

# Returns
- `rf`: (`::RF`) RF struct
"""
function get_RF(rfLibrary, shapeLibrary, Δt_rf, i)
    return get_RF(rfLibrary[i], shapeLibrary, Δt_rf)
end

pulseq_rf_adds_half_step(r::PulseqRFEvent) = r.time_shape_id <= 0

function get_RF(r::PulseqRFEvent, shapeLibrary, rf_raster_time)
    amplitude = r.amplitude
    mag_id = r.mag_id
    phase_id = r.phase_id
    time_shape_id = r.time_shape_id
    first_sample_offset = time_shape_id <= 0 ? rf_raster_time / 2 : 0.0
    delay = r.delay + first_sample_offset
    freq = r.freq
    phase = r.phase
    #Amplitude and phase waveforms
    if amplitude != 0 && mag_id != 0
        rfA = decompress_shape(shapeLibrary[mag_id]...)
        rfϕ = decompress_shape(shapeLibrary[phase_id]...)
        Nrf = shapeLibrary[mag_id][1] - 1
        rfAϕ = amplitude .* rfA .* exp.(1im * 2π * rfϕ)
    else
        rfAϕ = ComplexF64[0.0]
        Nrf = 0
    end
    #Creating timings
    if time_shape_id == 0 #no time waveform. Default time raster
        rfT = Nrf * rf_raster_time
    elseif time_shape_id == -1 #New in pulseq 1.5.x: no time waveform. 1/2 of the default time raster
        rfT = Nrf * rf_raster_time / 2
    else #time waveform
        rft = decompress_shape(shapeLibrary[time_shape_id]...)
        first_sample_offset = rft[1] * rf_raster_time
        delay += first_sample_offset # offset due to the shape starting at a non-zero value
        rfT = diff(rft) * rf_raster_time
    end
    center = isnothing(r.center) ? nothing : r.center - first_sample_offset
    return RF(rfAϕ, rfT, freq, delay; center, ϕ=phase, use=get_RF_use_from_char(Val(r.use)))
end

"""
    adc = get_ADC(adcLibrary, i)

Decode one Pulseq ADC event into a runtime `ADC`. Used by
[`decode_pulseq_libraries`](@ref).

# Arguments
- `adcLibrary`: (`::Dict{Int, PulseqADCEvent}`) the "adcLibrary" dictionary
- `i`: (`::Int64`) index of the adc in the block event

# Returns
- `adc`: (`::ADC`) ADC struct
"""
function get_ADC(adcLibrary, i)
    return get_ADC(adcLibrary[i])
end

function get_ADC(a::PulseqADCEvent)
    num = a.num
    dwell = a.dwell
    delay = a.delay + dwell / 2
    freq = a.freq
    phase = a.phase
    #Definition
    T = (num-1) * dwell
    return ADC(num,T,delay,freq,phase)
end

"""
    EXT = get_extension(extensionInstanceLibrary, extensionTypeLibrary, extensionSpecLibrary, i)

Reads the extension(s) for a block event in a Pulseq sequence file.

# Arguments
- `extensionInstanceLibrary`: (`::Dict{K, V}`) the extension library dictionary
- `extensionTypeLibrary`: (`::Dict{K, V}`) the extension type dictionary
- `extensionSpecLibrary`: (`::Dict{K, V}`) the extension specifications dictionary
- `i`: (`::Int64`) index of the extension in the block event

# Returns
- `EXT`: (`Vector{Extension}`) vector of Extension objects for the block event

# Details
The available extensions are currently contained in the file KomaMRIBase/src/datatypes/sequence/EXT.jl
"""
function get_extension(extensionInstanceLibrary, extensionTypeLibrary, extensionSpecLibrary, i)
    EXT = Extension[]
    haskey(extensionInstanceLibrary, i) || return EXT
    # type ref next_id
    # next_id of 0 terminates the list
    while i != 0
        entry = extensionInstanceLibrary[i]
        if haskey(extensionTypeLibrary, entry.type)
            push!(EXT, extensionSpecLibrary[entry.type][entry.ref])
        end
        i = entry.next_id
    end
    return EXT
end

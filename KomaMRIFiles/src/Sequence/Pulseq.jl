include("PulseqEvents.jl")

const QUANT_TOL = 1e-12
const DEFAULT_RASTER = PulseqRaster(Sequence(), Scanner())

"""
read_version Read the [VERSION] section of a sequence file.
   defs=read_version(fid) Read Pulseq version from file
   identifier of an open MR sequence file and return it
"""
function read_version(io)
    pulseq_version = VersionNumber(
        @scanf(readline(io), "major %i", Int)[end],
        @scanf(readline(io), "minor %i", Int)[end],
        @scanf(readline(io), "revision %i", Int)[end],
    )
    @assert pulseq_version.major == 1 "Unsupported version_major $(pulseq_version.major)"
    if     pulseq_version < v"1.2.0"
        @error "Unsupported Pulseq $(pulseq_version), only file format revision 1.2.0 and above are supported"
    elseif pulseq_version < v"1.3.1"
        @warn "Loading older Pulseq $(pulseq_version); some code may not function as expected"
    elseif pulseq_version >= v"1.5.0"
        @info "Pulseq $(pulseq_version) is supported, but Soft Delay, Rotation, and RF Shimming extensions are not yet included\n(see https://github.com/JuliaHealth/KomaMRI.jl/issues/714)" maxlog=1
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
    def = Dict{String,Any}()
    while true
        line = readline(io)
        line_split = String.(split(line))
        (length(line_split) > 0) || break #Break on white space
        key = line_split[1]
        value_string_array = line_split[2:end]
        parsed_array = [tryparse(Float64, s) === nothing ? s : tryparse(Float64, s) for s = value_string_array]
        def[key] = (length(parsed_array) == 1 && key != "RequiredExtensions") ? parsed_array[1] : parsed_array
    end
    if !haskey(def,"BlockDurationRaster")      def["BlockDurationRaster"]      = DEFAULT_RASTER.BlockDurationRaster end
    if !haskey(def,"GradientRasterTime")       def["GradientRasterTime"]       = DEFAULT_RASTER.GradientRasterTime end
    if !haskey(def,"RadiofrequencyRasterTime") def["RadiofrequencyRasterTime"] = DEFAULT_RASTER.RadiofrequencyRasterTime end
    if !haskey(def,"AdcRasterTime")            def["AdcRasterTime"]            = DEFAULT_RASTER.AdcRasterTime end
    return def
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
        line_split = String.(split(line))
        (length(line_split) > 0) || break #Break on white space
        firstchar = first(line_split[1])
        firstchar == '#' && continue
        key = line_split[1]
        value = join(line_split[2:end], " ")
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
    num_blocks = 0
    while true
        line = readline(io)
        isempty(line) && break
        blockEvents = [parse(Int, d) for d in eachsplit(line)]
        if length(blockEvents) != n_block_events
            error("Expected $n_block_events events but got $(length(blockEvents)).")
        end
        if blockEvents[1] != 0
            if pulseq_version <= v"1.2.1"
                # discard id and duration (duration treated below)
                append!(blocks, blockEvents[3:end])
                push!(blocks, 0) # append 0 for version compatibility (see assert below)
            else
                append!(blocks, blockEvents[3:end])
            end
            if pulseq_version >= v"1.4.0" # Explicit block duration (in units of blockDurationRaster)
                duration = blockEvents[2] * block_duration_raster
                push!(blockDurations, duration)
            else # Implicit block duration from delay ID
                push!(delayIDs_tmp, blockEvents[2])
            end
            num_blocks += 1
        end
    end
    reshaped_blocks = reshape(blocks, :, num_blocks)
    # we need 6 vals because we drop IDs and durations
    @assert size(reshaped_blocks, 1) == 6 "unexpected number of fields per block"
    return reshaped_blocks, blockDurations, delayIDs_tmp
end

apply_scale(scale::Number, value::Number) = scale * value
apply_scale(scale, value) = value

"""
    events = read_events(io, scale; format="%i "*"%f "^(length(scale)), event_library)

Read an event section of a sequence file.

# Arguments
- `io`: (`::IO`) input file stream
- `scale`: tuple of scale factors applied to parsed fields

# Keywords
- `format`: (`::String`) format string
- `event_library`: output dictionary for the parsed event library
- `constructor`: constructor used to normalize each parsed event row

# Returns
- `events`: parsed event library
"""
function read_events(io, scales; format="%i " * "%f " ^ length(scales), event_library, constructor=identity)
    event_length = length(scales) + 1
    fmt = Scanf.Format(format)
    args = Tuple(
        token == "%i" ? Int : token == "%f" ? Float64 : token == "%c" ? Char : String
        for token in split(format)
    )
    while true
        line = readline(io)
        isempty(line) && break
        r, data... = scanf(line, fmt, args...)
        r == event_length || break #Break if not all values read
        id = floor(Int, data[1])
        parsed = ntuple(i -> apply_scale(scales[i], data[i + 1]), length(scales))
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
        Tuple(KomaMRIBase.get_scale(ext_type));
        format="%i " * KomaMRIBase.get_scanf_format(ext_type),
        event_library=Dict{Int, Tuple}(),
    )
    extensionSpecLibrary[ext_id] = Dict(id => ext_type(values...) for (id, values) in raw_specs)
end
function read_extensions(io, ext_string, ext_type, ext_id, extensionTypeLibrary, extensionSpecLibrary, required_extensions)
    if ext_string in required_extensions
        error("Extension $ext_string is required by the sequence (RequiredExtensions: $required_extensions) but not supported by KomaMRI reader")
    else
        @warn "Ignoring unsupported extension: $ext_string"
        while true # Skip the extension specifications
            line = readline(io)
            isempty(line) && break
        end
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
        r, id = @scanf(line, "shape_id %i", Int)
        if r != 1 reset(io); break end # If the line is not a shape_id, reset the io and break the loop
        unmark(io) # Unmark the position after reading the line
        _, num_samples = @scanf(readline(io), "num_samples %i", Int)
        shape = Float64[]
        while true #Reading shape data
            data_point = tryparse(Float64, readline(io))
            !isnothing(data_point) || break #Break if no sample
            append!(shape, data_point)
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
compress_shape Compress a gradient or pulse shape.
   s=compress_shape(w) Compress the waveform using a run-length compression
   scheme on the derivative. This strategy encodes constant and linear
   waveforms with very few samples. Returns:
     num_samples - the number of samples in the uncompressed waveform
     data - containing the compressed waveform

   See also decompress_shape
"""
function compress_shape(w; forceCompression=false)
    num_samples = length(w)
    if !forceCompression && num_samples <= 4
        data = w[:]
    else
        quant_fac = 1e-7
        ws = w ./ quant_fac;
        datq = round.([ws[1]; diff(ws[:])])
        qerr = ws[:] .- cumsum(datq)
        qcor = [0; diff(round.(qerr))]
        datd = datq .+ qcor
        maskChanges = [true; diff(datd) .!= 0]
        vals = datd[maskChanges] .* quant_fac # Elements without repetitions
        k = findall([maskChanges; true]) # Indices of changes
        n = diff(k) # Number of repetitions
        # Encode in Pulseq format
        nExtra = convert(Vector{Any}, n .- 2)
        vals2 = convert(Vector{Any}, vals)
        vals2[nExtra .< 0] .= NaN;
        nExtra[nExtra .< 0] .= NaN;
        v = [vals'; vals2'; nExtra']
        v = convert(Vector{Float64}, v[isfinite.(v)])
        v[abs.(v) .<= 1e-10] .= 0
        # decide whether compression makes sense, otherwise store the original
        data = forceCompression || num_samples > length(v) ? v : w
    end
    return num_samples, data
end

"""
decompress_shape Decompress a gradient or pulse shape.
   w=decompress_shape(shape) Decompress the shape compressed with a run-length
   compression scheme on the derivative. Returns:
     num_samples - the number of samples in the uncompressed waveform
     data - containing the compressed waveform

   See also compress_shape
"""
function decompress_shape(num_samples, data; forceDecompression = false)
    dataPack = data
    dataPackLen = length(dataPack)
    numSamples = num_samples
    if !forceDecompression && numSamples == dataPackLen
        w = dataPack
    else
        w = zeros(numSamples) #pre-allocate the results
        #Decompression starts
        dataPackDiff = dataPack[2:end] .- dataPack[1:end-1]
        #when dataPackDiff == 0 the subsequent samples are equal ==> marker for
        #repeats (run-length encoding)
        dataPackMarkers = findall(dataPackDiff .== 0)
        countPack = 1       # counter 1: points to the current compressed sample
        countUnpack = 1     # counter 2: points to the current uncompressed sample
        for i = eachindex(dataPackMarkers)
            nextPack = dataPackMarkers[i] # careful, this index may have "false positives" , e.g. if the value 3 repeats 3 times, then we will have 3 3 3
            currUnpackSamples = nextPack - countPack
            if currUnpackSamples < 0 # this rejects false positives
                continue
            elseif currUnpackSamples > 0 #do we have an unpacked block to copy?
                w[countUnpack:(countUnpack+currUnpackSamples-1)] .= dataPack[countPack:(nextPack-1)]
                countPack += currUnpackSamples
                countUnpack += currUnpackSamples
            end
            # now comes the packed/repeated section
            rep = floor(Int, dataPack[countPack+2] + 2)
            w[countUnpack:(countUnpack+rep-1)] .= dataPack[countPack]
            countPack += 3
            countUnpack += rep
        end
        # samples left?
        if countPack <= dataPackLen
            # println("$dataPackLen $countPack $numSamples $countUnpack")
            @assert dataPackLen-countPack == numSamples-countUnpack "Unsuccessful unpacking of samples"
            # copy the rest of the shape, it is unpacked
            w[countUnpack:end] .= dataPack[countPack:end]
        end
        w = cumsum(w)
    end
    return w
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
    updated_grad_ids = BitSet()
    for (iB, eventIDs) in enumerate(eachcol(blockEvents))
        block_duration = blockDurations[iB]
        for iG in eachindex(grad_prev_last)
            g_id = eventIDs[1 + iG]
            g_id > 0 || continue
            update_first_last!(grad_prev_last, iG, eventLibraries.grad_library[g_id], block_duration, eventLibraries, updated_grad_ids, g_id)
        end
    end
end

function fix_legacy_trapezoids!(grad_library, grad_raster_time)
    for (i, grad) in grad_library
        grad isa TrapGradEvent || continue
        fix_rise = grad.rise == 0 && abs(grad.amplitude) == 0 && grad.flat > 0
        flat = fix_rise ? grad.flat - grad_raster_time : grad.flat
        fix_delay = grad.delay == 0 && abs(grad.amplitude) == 0 && flat > 0
        (fix_rise || fix_delay) || continue
        grad_library[i] = TrapGradEvent(
            grad.amplitude,
            fix_rise ? grad_raster_time : grad.rise,
            fix_delay ? flat - grad_raster_time : flat,
            grad.fall,
            fix_delay ? grad_raster_time : grad.delay,
        )
    end
end

function init_legacy_block_durations!(blockDurations, blockEvents, delayIDs_tmp, eventLibraries)
    resize!(blockDurations, size(blockEvents, 2))
    for (i, eventIDs) in enumerate(eachcol(blockEvents))
        delayID = delayIDs_tmp[i]
        delay = delayID > 0 ? eventLibraries.tmp_delay_library[delayID] : 0.0
        Gx, Gy, Gz, rf, add_half_Δt_rf, adc, _ = get_block(eventIDs, eventLibraries)
        blockDurations[i] = max(
            delay,
            dur(Gx),
            dur(Gy),
            dur(Gz),
            dur(rf) + add_half_Δt_rf * eventLibraries.radiofrequency_raster_time / 2,
            dur(adc),
        )
    end
end

function update_first_last!(grad_prev_last, iG, ::TrapGradEvent, block_duration, event_libraries, updated_grad_ids, g_id)
    grad_prev_last[iG] = 0.0
    return nothing
end

function update_first_last!(grad_prev_last, iG, g::ArbGradEvent, block_duration, event_libraries, updated_grad_ids, g_id)
    if g_id ∉ updated_grad_ids
        g.first = grad_prev_last[iG]
        g.last = 0.0
        waveform = g.amplitude * decompress_shape(event_libraries.shape_library[g.amp_shape_id]...)
        if g.time_shape_id != 0 # time-shaped case
            g.last = waveform[end]
        else # uniformly-shaped case
            odd_step1 = [g.first; 2 * waveform]
            odd_step2 = odd_step1 .* (mod.(1:length(odd_step1), 2) * 2 .- 1)
            waveform_odd_rest = cumsum(odd_step2) .* (mod.(1:length(odd_step2), 2) * 2 .- 1)
            g.last = waveform_odd_rest[end]
        end
        push!(updated_grad_ids, g_id)
    end
    grad_duration = dur(get_Grad(g, event_libraries.shape_library, event_libraries.gradient_raster_time))
    grad_prev_last[iG] = grad_duration + eps(Float64) < block_duration ? 0 : g.last
    return nothing
end

"""
    amp_shape, time_shape = simplify_waveforms(amp_shape, time_shape)

Simplifies the amplitude and time waveforms when time steps are uniform.
"""
@inline simplify_waveforms(amp_shape, time_shape::Real) = amp_shape, time_shape

function simplify_waveforms(amp_shape, time_shape::AbstractVector{<:Real})
    isempty(time_shape) && return amp_shape, time_shape
    all(==(time_shape[1]), time_shape) || return amp_shape, time_shape
    return amp_shape, sum(time_shape)
end

"""
    seq = get_seq_from_blocks(blockEvents, blockDurations, eventLibraries)

Sequence definition from several blocks, ordered by occurrence. Used internally by [`read_seq`](@ref).

# Arguments
- `blockEvents`: (`::Array{Int, 2}`) array of block event IDs
- `blockDurations`: (`::Vector{Float64}`) vector of block durations
- `eventLibraries`: (`::PulseqEventLibraries`) main dictionary of event libraries

# Returns
- `s`: (`::Sequence`) Sequence struct
"""
function get_seq_from_blocks(blockEvents, blockDurations, eventLibraries)
    num_blocks = size(blockEvents, 2)
    GR = Matrix{Grad}(undef, 3, num_blocks)
    RFs = Matrix{RF}(undef, 1, num_blocks)
    ADCs = Vector{ADC}(undef, num_blocks)
    EXTs = Vector{Vector{Extension}}(undef, num_blocks)

    for (i, eventIDs) in enumerate(eachcol(blockEvents))
        Gx, Gy, Gz, rf, _, adc, ext = get_block(eventIDs, eventLibraries)
        GR[1, i] = Gx
        GR[2, i] = Gy
        GR[3, i] = Gz
        RFs[1, i] = rf
        ADCs[i] = adc
        EXTs[i] = ext
    end

    # Definitions
    DEFs = Dict{String,Any}()

    # Sequence 
    return Sequence(GR, RFs, ADCs, blockDurations, EXTs, DEFs)
end

"""
    seq = read_seq(filename)

Returns the Sequence struct from a Pulseq file with `.seq` extension.

# Arguments
- `filename`: (`::String`) absolute or relative path of the sequence file `.seq`

# Returns
- `seq`: (`::Sequence`) Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_seq(seq)
```
"""
function read_seq(filename)
    @info "Loading sequence $(basename(filename)) ..."
    pulseq_version = v"0.0.0"
    gradLibrary = Dict{Int, GradEvent}()
    def = Dict{String, Any}()
    signature = nothing
    blockEvents = Array{Int, 2}(undef, 0, 0)
    blockDurations = Vector{Float64}(undef, 0)
    delayIDs_tmp = Int[]
    rfLibrary = Dict{Int, RFEvent}()
    adcLibrary = Dict{Int, ADCEvent}()
    tmp_delayLibrary = Dict{Int, Float64}()
    shapeLibrary = ShapeLibrary()
    extensionInstanceLibrary = Dict{Int, ExtensionInstanceEvent}()
    extensionTypeLibrary = Dict{Int, Type{<:Extension}}()
    extensionSpecLibrary = Dict{Int, Dict{Int, Extension}}()
    #Reading file and storing data
    open(filename) do io
        while !eof(io)
            section = readline(io)
            if isempty(section) || section[1] == '#'
                #skip useless line
            elseif section == "[DEFINITIONS]"
                def = read_definitions(io)
            elseif  section == "[VERSION]"
                pulseq_version = read_version(io)
            elseif  section == "[BLOCKS]"
                if pulseq_version == v"0.0.0"
                    @error "Pulseq file MUST include [VERSION] section prior to [BLOCKS] section"
                end
                block_duration_raster = get(def, "BlockDurationRaster", DEFAULT_RASTER.BlockDurationRaster)
                blockEvents, blockDurations, delayIDs_tmp = read_blocks(io, block_duration_raster, pulseq_version)
            elseif  section == "[RF]"
                if pulseq_version >= v"1.5.0"
                    rfLibrary = read_events(io, (1 / γ, 1, 1, 1, 1e-6, 1e-6, 1, 1, 1, 1, 1); format="%i " * "%f " ^ 10 * "%c ", event_library=rfLibrary, constructor=RFEvent) # this is 1.5.x format
                elseif pulseq_version >= v"1.4.0"
                    rfLibrary = read_events(io, (1 / γ, 1, 1, 1, 1e-6, 1, 1); event_library=rfLibrary, constructor=RFEvent) # this is 1.4.x format
                else
                    rfLibrary = read_events(io, (1 / γ, 1, 1, 1e-6, 1, 1); event_library=rfLibrary, constructor=RFEvent) # this is 1.3.x and below
                    # we will have to scan through the library later after all the shapes have been loaded
                end
            elseif  section == "[GRADIENTS]"
                if pulseq_version >= v"1.5.0"
                    gradLibrary = read_events(io, (1 / γ, 1 / γ, 1 / γ, 1, 1, 1e-6); event_library=gradLibrary, constructor=ArbGradEvent) # this is 1.5.x format
                elseif pulseq_version >= v"1.4.0"
                    gradLibrary = read_events(io, (1 / γ, 1, 1, 1e-6); event_library=gradLibrary, constructor=ArbGradEvent) # this is 1.4.x format
                else
                    gradLibrary = read_events(io, (1 / γ, 1, 1e-6); event_library=gradLibrary, constructor=ArbGradEvent) # this is 1.3.x and below
                end
            elseif  section == "[TRAP]"
                gradLibrary = read_events(io, (1 / γ, 1e-6, 1e-6, 1e-6, 1e-6); event_library=gradLibrary, constructor=TrapGradEvent)
            elseif  section == "[ADC]"
                if pulseq_version >= v"1.5.0"
                    adcLibrary = read_events(io, (1, 1e-9, 1e-6, 1, 1, 1, 1, 1); event_library=adcLibrary, constructor=ADCEvent) # this is 1.5.x format
                else
                    adcLibrary = read_events(io, (1, 1e-9, 1e-6, 1, 1); event_library=adcLibrary, constructor=ADCEvent) # this is 1.4.x and below
                end
            elseif  section == "[DELAYS]"
                if pulseq_version >= v"1.4.0"
                    @error "Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section"
                end
                tmp_delayLibrary = read_events(io, (1e-6,); event_library=tmp_delayLibrary, constructor=first);
            elseif  section == "[SHAPES]"
                shapeLibrary = read_shapes(io, (pulseq_version.major == 1 && pulseq_version.minor < 4))
            elseif  section == "[EXTENSIONS]"
                extensionInstanceLibrary = read_events(io, (1, 1, 1); format="%i " * "%i " ^ 3, event_library=extensionInstanceLibrary, constructor=ExtensionInstanceEvent)
            elseif  section == "[SIGNATURE]"
                signature = read_signature(io)
            elseif startswith(section, "extension")
                ext = section[11:end]
                ext_string = split(ext, " ")[1]
                ext_type   = KomaMRIBase.get_EXT_type_from_symbol(Val(Symbol(ext_string)))
                ext_id     = parse(Int, split(ext, " ")[2])
                if !haskey(def,"RequiredExtensions") def["RequiredExtensions"] = String[]  end
                read_extensions(io, ext_string, ext_type, ext_id, extensionTypeLibrary, extensionSpecLibrary, def["RequiredExtensions"])
            else
                @error "Unknown section code: $section"
            end
        end
    end

    # fix trapezoidal gradients imported from older versions
    if pulseq_version < v"1.4.0"
        grad_raster_time = get(def, "gradRasterTime", get(def, "GradientRasterTime", DEFAULT_RASTER.GradientRasterTime))
        fix_legacy_trapezoids!(gradLibrary, grad_raster_time)
    end
    verify_signature!(filename, signature; pulseq_version=pulseq_version)
    isempty(def) && (def = merge_definitions_with_raster!(Dict{String, Any}(), DEFAULT_RASTER))
    #Sequence libraries (basically everything except the blocks)
    eventLibraries = PulseqEventLibraries(
        gradLibrary,
        rfLibrary,
        adcLibrary,
        tmp_delayLibrary,
        shapeLibrary,
        extensionInstanceLibrary,
        extensionTypeLibrary,
        extensionSpecLibrary,
        def,
    )

    if pulseq_version < v"1.4.0"
        init_legacy_block_durations!(blockDurations, blockEvents, delayIDs_tmp, eventLibraries)
    end
    # Add first and last points for gradients #320 for version <= 1.4.2
    if pulseq_version < v"1.5.0"
        fix_first_last_grads!(blockEvents, blockDurations, eventLibraries)
    end

    seq = get_seq_from_blocks(blockEvents, blockDurations, eventLibraries)

    # Final details
    #Temporary hack
    seq.DEF = merge(eventLibraries.definitions, seq.DEF)
    # Koma specific details for reconstrucion
    seq.DEF["FileName"] = basename(filename)
    seq.DEF["PulseqVersion"] = pulseq_version
    seq.DEF["signature"] = signature
    # Guessing recon dimensions
    seq.DEF["Nx"] = trunc(Int64, get(seq.DEF, "Nx", maximum(adc.N for adc = seq.ADC)))
    seq.DEF["Nz"] = trunc(Int64, get(seq.DEF, "Nz", length(unique(seq.RF.Δf))))
    seq.DEF["Ny"] = trunc(Int64, get(seq.DEF, "Ny", sum(map(is_ADC_on, seq)) ÷ seq.DEF["Nz"]))
    #Koma sequence
    return seq
end

#To Sequence
"""
    grad = get_Grad(gradLibrary, shapeLibrary, Δt_gr, i)

Reads the gradient. It is used internally by [`get_block`](@ref).

# Arguments
- `gradLibrary`: (`::Dict{Int, GradEvent}`) the "gradLibrary" dictionary
- `shapeLibrary`: (`::ShapeLibrary`) the "shapeLibrary" dictionary
- `Δt_gr`: (`::Float64`, `[s]`) gradient raster time
- `i`: (`::Int64`) index of the axis in the block event

# Returns
- `grad`: (::Grad) Gradient struct
"""
function get_Grad(gradLibrary, shapeLibrary, Δt_gr, i)
    haskey(gradLibrary, i) || return Grad(0.0, 0.0)
    return get_Grad(gradLibrary[i], shapeLibrary, Δt_gr)
end

get_Grad(grad::TrapGradEvent, shapeLibrary, Δt_gr) = Grad(grad.amplitude, grad.flat, grad.rise, grad.fall, grad.delay, 0.0, 0.0)

function get_Grad(grad::ArbGradEvent, shapeLibrary, Δt_gr)
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
    gA, gT = simplify_waveforms(gA, gT)
    return Grad(gA, gT, rise, fall, delay, first_grads, last_grads)
end

"""
    rf, add_half_Δt_rf = get_RF(rfLibrary, shapeLibrary, Δt_rf, i)

Reads the RF. It is used internally by [`get_block`](@ref).

# Arguments
- `rfLibrary`: (`::Dict{Int, RFEvent}`) the "rfLibrary" dictionary
- `shapeLibrary`: (`::ShapeLibrary`) the "shapeLibrary" dictionary
- `Δt_rf`: (`::Float64`, `[s]`) RF raster time
- `i`: (`::Int64`) index of the RF in the block event

# Returns
- `rf`: (`::RF`) RF struct
- `add_half_Δt_rf`: (`::Bool`) whether the block extent needs the trailing half-raster correction
"""
function get_RF(rfLibrary, shapeLibrary, Δt_rf, i)
    rf = RF(0.0,0.0)
    haskey(rfLibrary, i) || return rf, false
    #Unpacking
    r = rfLibrary[i]
    amplitude = r.amplitude
    mag_id = r.mag_id
    phase_id = r.phase_id
    time_shape_id = r.time_shape_id
    add_half_Δt_rf = time_shape_id <= 0
    delay = r.delay + add_half_Δt_rf * Δt_rf / 2
    freq = r.freq
    phase = r.phase
    #Amplitude and phase waveforms
    if amplitude != 0 && mag_id != 0
        rfA = decompress_shape(shapeLibrary[mag_id]...)
        rfϕ = decompress_shape(shapeLibrary[phase_id]...)
        Nrf = shapeLibrary[mag_id][1] - 1
        rfAϕ = amplitude .* rfA .* exp.(1im*(2π*rfϕ .+ phase))
    else
        rfAϕ = ComplexF64[0.0]
        Nrf = 0
    end
    #Creating timings
    if time_shape_id == 0 #no time waveform. Default time raster
        rfT = Nrf * Δt_rf
    elseif time_shape_id == -1 #New in pulseq 1.5.x: no time waveform. 1/2 of the default time raster
        rfT =  Nrf * Δt_rf / 2
    else #time waveform
        rft = decompress_shape(shapeLibrary[time_shape_id]...)
        delay += rft[1] * Δt_rf # offset due to the shape starting at a non-zero value
        rfT = diff(rft) * Δt_rf
    end
    rfAϕ, rfT = simplify_waveforms(rfAϕ, rfT)
    if isnothing(r.center)
        rf = RF(rfAϕ, rfT, freq, delay)
    else
        rf = RF(rfAϕ, rfT, freq, delay, r.center, KomaMRIBase.get_RF_use_from_char(Val(r.use)))
    end
    return rf, add_half_Δt_rf
end

"""
    adc = get_ADC(adcLibrary, i)

Reads the ADC. It is used internally by [`get_block`](@ref).

# Arguments
- `adcLibrary`: (`::Dict{Int, ADCEvent}`) the "adcLibrary" dictionary
- `i`: (`::Int64`) index of the adc in the block event

# Returns
- `adc`: (`::ADC`) ADC struct
"""
function get_ADC(adcLibrary, i)
    adc = ADC(0, 0)
    haskey(adcLibrary, i) || return adc
    #Unpacking
    a = adcLibrary[i]
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
    Gx, Gy, Gz, rf, add_half_Δt_rf, adc, ext = get_block(eventIDs, eventLibraries)

Decode one Pulseq block into gradient, RF, ADC, and extension events.

# Arguments
- `eventIDs`: (`::Vector{Int}`) event IDs for one block
- `eventLibraries`: (`::PulseqEventLibraries`) main dictionary of event libraries

# Returns
- `Gx`, `Gy`, `Gz`: (`::Grad`) decoded gradient events
- `rf`: (`::RF`) decoded RF event
- `add_half_Δt_rf`: (`::Bool`) whether the block extent needs the trailing half-raster correction
- `adc`: (`::ADC`) decoded ADC event
- `ext`: (`::Vector{Extension}`) decoded extension events
"""
function get_block(eventIDs, eventLibraries)
    irf, igx, igy, igz, iadc, iext = eventIDs
    Gx = get_Grad(eventLibraries.grad_library, eventLibraries.shape_library, eventLibraries.gradient_raster_time, igx)
    Gy = get_Grad(eventLibraries.grad_library, eventLibraries.shape_library, eventLibraries.gradient_raster_time, igy)
    Gz = get_Grad(eventLibraries.grad_library, eventLibraries.shape_library, eventLibraries.gradient_raster_time, igz)
    rf, add_half_Δt_rf = get_RF(eventLibraries.rf_library, eventLibraries.shape_library, eventLibraries.radiofrequency_raster_time, irf)
    adc = get_ADC(eventLibraries.adc_library, iadc)
    ext = get_extension(eventLibraries.extension_instance_library, eventLibraries.extension_type_library, eventLibraries.extension_spec_library, iext)
    return Gx, Gy, Gz, rf, add_half_Δt_rf, adc, ext
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
            push!(EXT, deepcopy(extensionSpecLibrary[entry.type][entry.ref]))
        end
        i = entry.next_id
    end
    return EXT
end

# ----------------- Signature-related functions ------------------
function supported_signature_digest(algorithm::AbstractString, payload::Vector{UInt8})
    alg = lowercase(strip(algorithm))
    digest = if alg == "md5"
        md5(payload)
    elseif alg == "sha1"
        sha1(payload)
    elseif alg in ("sha2", "sha256")
        sha256(payload)
    else
        throw(ArgumentError("Unsupported signature algorithm '$algorithm'. Supported algorithms: md5, sha1, sha256."))
    end
    return lowercase(bytes2hex(digest))
end

function parse_signature_section(section::AbstractString)
    io = IOBuffer(section)
    readline(io) # consume "[SIGNATURE]" header
    signature_type = nothing
    signature_hash = nothing
    while !eof(io)
        line = strip(readline(io))
        isempty(line) && continue
        startswith(line, '#') && continue
        parts = split(line)
        isempty(parts) && continue
        key = parts[1]
        value = join(parts[2:end], " ")
        if key == "Type"
            signature_type = strip(value)
        elseif key == "Hash"
            signature_hash = lowercase(replace(strip(value), " " => ""))
        end
    end
    return isnothing(signature_type) && isnothing(signature_hash) ? nothing : (type = signature_type, hash = signature_hash)
end

function verify_signature!(filename::String, signature::Nothing; pulseq_version::VersionNumber=v"1.4.0")
    @warn "Pulseq [SIGNATURE] section is missing; skipping verification."
    return nothing
end

function verify_signature!(filename::String, signature::NamedTuple; pulseq_version::VersionNumber=v"1.4.0")
    # Don't error, just warn - the file can still be used without the signature or with it being incorrect
    sig_type = signature.type
    sig_hash = signature.hash
    # Read file as bytes to avoid any line ending normalization
    file_bytes = read(filename)
    # Find [SIGNATURE] in bytes
    sig_marker = b"[SIGNATURE]"
    sig_pos = findfirst(sig_marker, file_bytes)
    isnothing(sig_pos) && begin
        @warn "Signature section expected but not found when verifying Pulseq file." filename
        return
    end
    sig_start = first(sig_pos)
    # Extract payload based on version
    # For Pulseq < 1.4.0 (e.g., JEMRIS), the newline before [SIGNATURE] is part of the payload
    # For Pulseq >= 1.4.0, the newline before [SIGNATURE] is part of the signature and should be excluded
    # However, different implementations may handle this differently, so we try multiple approaches
    payload_end = sig_start - 1
    expected_hash = lowercase(replace(sig_hash, " " => ""))
    if pulseq_version < v"1.4.0"
        # Include the newline before [SIGNATURE] in the payload (JEMRIS format)
        payload_bytes = payload_end > 0 ? file_bytes[1:payload_end] : UInt8[]
        computed_hash = supported_signature_digest(sig_type, payload_bytes)
    else
        # For version >= 1.4.0, try multiple approaches to handle different implementations
        # 1. Exclude only the last newline (spec-compliant)
        if payload_end > 0 && file_bytes[payload_end] in (UInt8('\n'), UInt8('\r'))
            payload_bytes_excluding = file_bytes[1:(payload_end - 1)]
        else
            payload_bytes_excluding = payload_end > 0 ? file_bytes[1:payload_end] : UInt8[]
        end
        # 2. Include the last newline (some implementations like MATLAB)
        payload_bytes_including = file_bytes[1:payload_end]
        # 3. Exclude all consecutive newlines before [SIGNATURE] (some implementations)
        last_non_nl = payload_end
        while last_non_nl > 0 && file_bytes[last_non_nl] in (UInt8('\n'), UInt8('\r'))
            last_non_nl -= 1
        end
        payload_bytes_no_nl = last_non_nl > 0 ? file_bytes[1:last_non_nl] : UInt8[]
        # Try all approaches and use the one that matches
        computed_hash_excluding = supported_signature_digest(sig_type, payload_bytes_excluding)
        computed_hash_including = supported_signature_digest(sig_type, payload_bytes_including)
        computed_hash_no_nl = supported_signature_digest(sig_type, payload_bytes_no_nl)
        if computed_hash_excluding == expected_hash
            computed_hash = computed_hash_excluding
            payload_bytes = payload_bytes_excluding
        elseif computed_hash_including == expected_hash
            computed_hash = computed_hash_including
            payload_bytes = payload_bytes_including
        elseif computed_hash_no_nl == expected_hash
            computed_hash = computed_hash_no_nl
            payload_bytes = payload_bytes_no_nl
        else
            # None matches, use the spec-compliant one for error message
            computed_hash = computed_hash_excluding
            payload_bytes = payload_bytes_excluding
        end
    end
    # Parse signature section for comparison (read as string for parsing)
    file_text = String(file_bytes)
    sig_section_start = sig_start
    sig_section = file_text[sig_section_start:end]
    parsed = parse_signature_section(sig_section)
    isnothing(parsed) && begin
        @warn "Failed to parse Pulseq signature section; skipping verification." filename
        return
    end
    parsed_type = parsed.type
    parsed_hash = parsed.hash
    if !isnothing(parsed_type) && !isempty(parsed_type) && lowercase(parsed_type) != lowercase(sig_type)
        @warn "Pulseq signature Type mismatch between parser and metadata. Using '$sig_type'." filename parsed_type
    end
    if !isnothing(parsed_hash) && !isempty(parsed_hash)
        parsed_hash_norm = lowercase(replace(parsed_hash, " " => ""))
        if parsed_hash_norm != lowercase(replace(sig_hash, " " => ""))
            @warn "Pulseq signature Hash mismatch between parser and metadata. Using metadata value." filename
        end
    end
    if computed_hash != expected_hash
        @warn "Pulseq signature verification failed for $(basename(filename)). Expected $(expected_hash), computed $(computed_hash). The file may have been modified or generated with a different implementation." filename
    end
end

# ----------------- Write Pulseq -----------------

"""
    collect_pulseq_assets(seq::Sequence, raster::PulseqRaster) -> (blocks, event_libraries)

Create the Pulseq export dictionaries required to serialize `seq` into the Pulseq file format.
This function is responsible for deduplicating reusable objects (RF, gradients, shapes, etc.)
and for translating each sequence block into integer lookups expected by the specification.
"""
function collect_pulseq_assets(seq::Sequence, raster::PulseqRaster)
    blocks = Dict{Int,PulseqBlock}()
    rf_library = Dict{Int,RFEvent}()
    grad_library = Dict{Int,GradEvent}()
    adc_library = Dict{Int,ADCEvent}()
    tmp_delay_library = Dict{Int,Float64}()
    shape_library = ShapeLibrary()
    extension_instance_library = Dict{Int,ExtensionInstanceEvent}()
    extension_type_library = Dict{Int,Type{<:Extension}}()
    extension_spec_library = Dict{Int,Dict{Int,Extension}}()
    definitions = Dict{String, Any}(seq.DEF)
    
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
            ext_id = register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext_vec)
            push!(extension_vectors, ext_vec)
            push!(extension_vector_ids, ext_id)
        end
    end

    # Second step: process all blocks and use pre-registered extension IDs.
    for (idx, block) in enumerate(seq)
        rf = block.RF[1]
        rf_id = register_rf!(rf_library, shape_library, rf.A, rf.T, rf.Δf, rf.delay, rf.center, rf.use, raster)

        gx, gy, gz = block.GR
        gx_id = register_grad!(grad_library, shape_library, gx.A, gx.T, gx.rise, gx.fall, gx.delay, gx.first, gx.last, raster)
        gy_id = register_grad!(grad_library, shape_library, gy.A, gy.T, gy.rise, gy.fall, gy.delay, gy.first, gy.last, raster)
        gz_id = register_grad!(grad_library, shape_library, gz.A, gz.T, gz.rise, gz.fall, gz.delay, gz.first, gz.last, raster)

        adc = block.ADC[1]
        adc_id, adc_dur = register_adc!(adc_library, adc.N, adc.T, adc.delay, adc.Δf, adc.ϕ)

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
        blocks[idx] = PulseqBlock(idx, duration, rf_id, gx_id, gy_id, gz_id, adc_id, ext_id)
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
function register_rf!(rf_library, shape_library, A, T, Δf, delay, center, use, raster::PulseqRaster)
    iszero(maximum(abs.(A))) && return 0
    Δt_rf = raster.RadiofrequencyRasterTime
    mag_id, phase_id, time_id = register_rf_shapes!(shape_library, A, T, Δt_rf; compress=true)
    amp = γ * maximum(abs.(A)) # from T to Hz (nucleus-dependent)
    delay = delay - (time_id == 0) * Δt_rf / 2
    phase  = mod.(angle.(A), 2π)[1]
    use_char = KomaMRIBase.get_char_from_RF_use(use)
    rf_event = RFEvent(amp, mag_id, phase_id, time_id, center, delay, 0.0, 0.0, Δf, phase, use_char)
    return _store_event!(rf_library, rf_event)
end

# A and T are numbers (pulse waveform)
function register_rf_shapes!(shapes, A, T, Δrf; compress = true)
    mag_id   = _store_shape!(shapes, [1.0, 1.0]; compress=compress)
    phase_id = _store_shape!(shapes, [0.0, 0.0]; compress=compress, phase_shape=true)
    time_id  = _store_shape!(shapes, [0.0, T] ./ Δrf; compress=compress)
    return mag_id, phase_id, time_id
end
# A is a vector (uniformly-sampled waveform)
function register_rf_shapes!(shapes, A::Vector, T, Δrf; compress = true)
    n_samples   = length(A)
    mag_id      = _store_shape!(shapes, abs.(A) ./ maximum(abs.(A)); compress=compress)
    phase_shape = mod.(angle.(A), 2π) / 2π
    phase_id    = _store_shape!(shapes, phase_shape .- phase_shape[1]; compress=compress, phase_shape=true)
    t_vector    = collect(range(0, T, length=n_samples)) ./ Δrf
    time_id     = _store_shape!(shapes, t_vector; compress=compress)
    return mag_id, phase_id, time_id
end
# A and T are vectors (time-shaped waveform)
function register_rf_shapes!(shapes, A::Vector, T::Vector, Δrf; compress = true)
    mag_id      = _store_shape!(shapes, abs.(A) ./ maximum(abs.(A)); compress=compress)
    phase_shape = mod.(angle.(A), 2π) / 2π
    phase_id    = _store_shape!(shapes, phase_shape .- phase_shape[1]; compress=compress, phase_shape=true)
    t_vector    = cumsum(vcat(0.0, T ./ Δrf))
    time_id     = _store_shape!(shapes, t_vector; compress=compress)
    return mag_id, phase_id, time_id
end

"""
    id = register_grad!(grad_library, shape_library, A, T, rise, fall, delay, first, last, raster)
"""
# Time-shaped waveform
function register_grad!(grad_library, shape_library, A::Vector, T::Vector, rise, fall, delay, first, last, raster::PulseqRaster)
    iszero(maximum(abs.(A))) && return 0
    if (iszero(rise) && iszero(fall))
        shape_id = _store_shape!(shape_library, A ./ maximum(abs.(A)); compress=true)
        t_vector = cumsum([0; T]) ./ raster.GradientRasterTime
        time_id  = _store_shape!(shape_library, t_vector; compress=true)
        amp   = γ * maximum(abs.(A)) # from T/m to Hz/m (nucleus-dependent)
        first = γ * first # from T/m to Hz/m 
        last  = γ * last  # from T/m to Hz/m
        return _store_event!(grad_library, ArbGradEvent(amp, first, last, shape_id, time_id, delay))
    end
    return register_grad!(grad_library, shape_library, [first; A; last], [rise; T; fall], 0, 0, delay, first, last, raster)
end
# Uniformly-sampled waveform
function register_grad!(grad_library, shape_library, A::Vector, T::Number, rise, fall, delay, first, last, raster::PulseqRaster)
    iszero(maximum(abs.(A))) && return 0
    Δgr = raster.GradientRasterTime
    intervals = diff(collect(range(0, T, length=length(A))))
    if (rise == Δgr/2 && fall == Δgr/2) & all(intervals .≈ Δgr)
        shape_id  = _store_shape!(shape_library, A ./ maximum(abs.(A)); compress=true)
        amp       = γ * maximum(abs.(A)) # from T/m to Hz/m (nucleus-dependent)
        first     = γ * first # from T/m to Hz/m 
        last      = γ * last  # from T/m to Hz/m
        return _store_event!(grad_library, ArbGradEvent(amp, first, last, shape_id, 0, delay))
    end
    return register_grad!(grad_library, shape_library, [first; A; last], [rise; intervals; fall], 0, 0, delay, first, last, raster)
end
# Trapezoidal waveform
function register_grad!(grad_library, shape_library, A::Number, T::Number, rise, fall, delay, first, last, raster::PulseqRaster)
    iszero(A) && return 0
    if (iszero(first) && iszero(last)) 
        return _store_event!(grad_library, TrapGradEvent(γ * A, rise, T, fall, delay))
    end
    return register_grad!(grad_library, shape_library, [first, A, A, last], [rise, T, fall], 0, 0, delay, first, last, raster)
end 

"""
    id = register_adc!(adc_library, N, T, delay, Δf, ϕ)

In Pulseq the ADC event starts at a dwell-time edge, but the first sample is taken at dwell/2
(see [Pulseq time and shape specification](https://pulseq.github.io/pulseq_shapes_and_times.pdf)).
Koma uses delay = time to first sample, so we store pulseq_delay = delay - dwell_s/2. For exact
round-trip (write → read) the Koma ADC delay must be ≥ dwell_s/2; otherwise it is clamped and a
warning is emitted.
"""
function register_adc!(adc_library, N, T, delay, Δf, ϕ)
    iszero(N) && return (0, 0.0)
    dwell_s = N == 1 ? T : T / (N - 1)
    delay_s = delay - dwell_s / 2
    adc_event = ADCEvent(N, dwell_s, delay_s, 0.0, 0.0, Δf, ϕ, 0)
    return _store_event!(adc_library, adc_event), delay_s + T + dwell_s
end

"""
    id = register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext)
"""
function register_ext!(extension_instance_library, extension_type_library, extension_spec_library, ext::Vector{Extension})
    length(ext) == 0 && return 0
    instance_ids = Int[]
    for e in ext
        ext_id = _store_event!(extension_type_library, typeof(e))
        if !haskey(extension_spec_library, ext_id)
            extension_spec_library[ext_id] = Dict{Int,Extension}()
        end
        ref = _store_event!(extension_spec_library[ext_id], e)
        instance_id = _store_event!(extension_instance_library, ExtensionInstanceEvent((ext_id, ref, 0)))
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

function _store_shape!(
    shapes::Dict{Int,Tuple{Int,Vector{Float64}}},
    samples::Vector{Float64};
    compress::Bool = true,
    phase_shape::Bool = false
)
    payload = compress ? compress_shape(samples) : (length(samples), samples)
    for (k, existing) in shapes
        if phase_shape
            existing_phase = existing[2] .- existing[2][1]
            if safe_equal_angles(existing_phase, payload[2]) && existing[1] == payload[1] return k end
        else
            if existing == payload return k end
        end
    end
    new_id = maximum(keys(shapes); init=0) + 1
    shapes[new_id] = payload
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

function emit_blocks_section!(io::IO, blocks)
    write(io, "# Format of blocks:\n")
    write(io, "# NUM DUR RF  GX  GY  GZ  ADC  EXT\n")
    write(io, "[BLOCKS]\n")
    isempty(blocks) && return
    block_ids = sort(collect(keys(blocks)))
    for id in block_ids
        block = blocks[id]
        values = [id, block.duration, block.rf, block.gx, block.gy, block.gz, block.adc, block.ext]
        max_lengths = zeros(Int, length(values))
        for (i, val) in enumerate(values)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, values, max_lengths)
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
    rf_ids = sort(collect(keys(event_libraries.rf_library)))
    for id in rf_ids
        rf_data = event_libraries.rf_library[id]
        center_us = isnothing(rf_data.center) ? 0.0 : rf_data.center * 1e6
        delay_us = round(Int, rf_data.delay * 1e6)
        values = [id, rf_data.amplitude, rf_data.mag_id, rf_data.phase_id, rf_data.time_shape_id, center_us, delay_us, rf_data.freq_ppm, rf_data.phase_ppm, rf_data.freq, rf_data.phase, rf_data.use]
        max_lengths = zeros(Int, 12)
        for (i, val) in enumerate(values)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, values, max_lengths)
    end
end

function emit_gradients_section!(io::IO, event_libraries)
    write(io, "\n# Format of arbitrary gradient events:\n")
    write(io, "# id      amp      first      last  shape_id  time_id  delay\n")
    write(io, "# ..     Hz/m       Hz/m      Hz/m        ..       ..     us\n")
    write(io, "[GRADIENTS]\n")
    !any(g -> g isa ArbGradEvent, values(event_libraries.grad_library)) && return
    grad_ids = sort(collect(keys(event_libraries.grad_library)))
    for id in grad_ids
        grad_data = event_libraries.grad_library[id]
        if grad_data isa ArbGradEvent
            vals = (id, grad_data.amplitude, grad_data.first, grad_data.last, grad_data.amp_shape_id, grad_data.time_shape_id, round(Int, grad_data.delay * 1e6))
            max_lengths = zeros(Int, 7)
            for (i, val) in enumerate(vals)
                str = _format_value(val)
                max_lengths[i] = max(max_lengths[i], length(str))
            end
            _format_row(io, vals, max_lengths)
        end
    end
end

function emit_trap_section!(io::IO, event_libraries)
    write(io, "\n# Format of trapezoid gradient events:\n")
    write(io, "# id      amp      rise  flat  fall  delay\n")
    write(io, "# ..     Hz/m        us    us    us     us\n")
    write(io, "[TRAP]\n")
    !any(g -> g isa TrapGradEvent, values(event_libraries.grad_library)) && return
    trap_ids = sort(collect(keys(event_libraries.grad_library)))
    for id in trap_ids
        trap_data = event_libraries.grad_library[id]
        if trap_data isa TrapGradEvent
            vals = (id, trap_data.amplitude, round(Int, trap_data.rise * 1e6), round(Int, trap_data.flat * 1e6), round(Int, trap_data.fall * 1e6), round(Int, trap_data.delay * 1e6))
            max_lengths = zeros(Int, 6)
            for (i, val) in enumerate(vals)
                str = _format_value(val)
                max_lengths[i] = max(max_lengths[i], length(str))
            end
            _format_row(io, vals, max_lengths)
        end
    end
end

function emit_adc_section!(io::IO, event_libraries)
    write(io, "\n# Format of ADC events:\n")
    write(io, "# id  num  dwell  delay  freq_ppm  phase_ppm  freq  phase  phase_id\n")
    write(io, "# ..   ..     ns     us       ppm        ppm    Hz    rad        ..\n")
    write(io, "[ADC]\n")
    isempty(event_libraries.adc_library) && return
    adc_ids = sort(collect(keys(event_libraries.adc_library)))
    for id in adc_ids
        adc_data = event_libraries.adc_library[id]
        vals = (id, adc_data.num, round(Int, adc_data.dwell * 1e9), round(Int, adc_data.delay * 1e6), adc_data.freq_ppm, adc_data.phase_ppm, adc_data.freq, adc_data.phase, adc_data.phase_id)
        max_lengths = zeros(Int, 9)
        for (i, val) in enumerate(vals)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, vals, max_lengths)
    end
end

function emit_extension_section!(io::IO, event_libraries)
    write(io, "\n# Format of extension events:\n")
    write(io, "# id  type  ref  next_id\n")
    write(io, "# Extension list is followed by extension specifications\n")
    write(io, "[EXTENSIONS]\n")
    isempty(event_libraries.extension_instance_library) && return
    ext_ids = sort(collect(keys(event_libraries.extension_instance_library)))
    for id in ext_ids
        ext_data = event_libraries.extension_instance_library[id]
        vals = (id, ext_data.type, ext_data.ref, ext_data.next_id)
        max_lengths = length.(digits.(vals))
        _format_row(io, vals, max_lengths)
    end
    write(io, "\n")
    type_ids = sort(collect(keys(event_libraries.extension_type_library)))
    for id in type_ids
        ext_type = event_libraries.extension_type_library[id]
        write(io, KomaMRIBase.extension_type_header(ext_type))
        write(io, "extension $(string(KomaMRIBase.get_symbol_from_EXT_type(ext_type))) $id\n")
        spec_ids = sort(collect(keys(event_libraries.extension_spec_library[id])))
        for spec_id in spec_ids
            v = event_libraries.extension_spec_library[id][spec_id]
            field_values = Tuple(getfield(v, f) for f in fieldnames(typeof(v)))
            scaled_values = [d isa Number ? s*d : d for (s, d) in zip(KomaMRIBase.get_scale(typeof(v)), field_values)]   
            vals = (spec_id, scaled_values...)
            max_lengths = zeros(Int, length(vals))
            for (i, val) in enumerate(vals)
                str = _format_value(val)
                max_lengths[i] = max(max_lengths[i], length(str))
            end
            _format_row(io, vals, max_lengths)
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

# Helper function to format a value as string based on format type
_format_value(val) = val isa Float64 && mod(val, 1.0) != 0.0 ? @sprintf("%.8f", val) : string(val)

# Helper function to format a table row with automatic column alignment (right-aligned)
function _format_row(io, values, max_lengths)
    for (i, (val, max_len)) in enumerate(zip(values, max_lengths))
        str = _format_value(val)
        # Right-align: add spaces before the value to pad to max_len
        num_spaces = max(0, max_len - length(str))
        write(io, repeat(" ", num_spaces))
        write(io, str)
        # Add spacing between columns (except after last column)
        if i < length(values)
            write(io, "   ")  # Base spacing between columns
        end
    end
    write(io, "\n")
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
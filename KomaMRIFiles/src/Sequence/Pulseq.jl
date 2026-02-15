const DEFAULT_DEFINITIONS = Dict("BlockDurationRaster"      => 1e-5, 
                                 "GradientRasterTime"       => 1e-5, 
                                 "RadiofrequencyRasterTime" => 1e-6, 
                                 "AdcRasterTime"            => 1e-7)

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
        def[key] = length(parsed_array) == 1 ? parsed_array[1] : parsed_array
    end
    #Default values
    if !haskey(def,"BlockDurationRaster")       def["BlockDurationRaster"] = DEFAULT_DEFINITIONS["BlockDurationRaster"] end
    if !haskey(def,"GradientRasterTime")        def["GradientRasterTime"] = DEFAULT_DEFINITIONS["GradientRasterTime"] end
    if !haskey(def,"RadiofrequencyRasterTime")  def["RadiofrequencyRasterTime"] = DEFAULT_DEFINITIONS["RadiofrequencyRasterTime"] end
    if !haskey(def,"AdcRasterTime")             def["AdcRasterTime"] = DEFAULT_DEFINITIONS["AdcRasterTime"] end
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
   library=read_blocks(fid) Read blocks from file identifier of an
   open MR sequence file and return the event table.
"""
function read_blocks(io, pulseq_version)
    eventTable = Dict{Int64, Vector{Int64}}()
    while true
        NumberBlockEvents = pulseq_version <= v"1.2.1" ? 7 : 8
        read_event = readline(io)
        !isempty(read_event) || break
        blockEvents = parse.(Int64, split(read_event))
        if blockEvents[1] != 0
            if pulseq_version <= v"1.2.1"
                eventTable[blockEvents[1]] = Int64[blockEvents[2:end]...; 0]
            else
                eventTable[blockEvents[1]] = Int64[blockEvents[2:end]...]
            end
        end
        length(blockEvents) == NumberBlockEvents || break #Break on white space
    end
    return eventTable
end

"""
    events = read_events(io, scale; type=nothing, format="%i "*"%f "^(length(scale)), eventLibrary=Dict())

Read an event section of a sequence file.

# Arguments
- `io`: (`::IO`) input file stream
- `scale`: (`::Vector{Float64}`) scale vector

# Keywords
- `type`: (`::Union{Nothing, Char}`) type of the event. Only used for gradients ('t' for trapezoidal and 'g' for arbitrary)
- `format`: (`::String`) format string
- `eventLibrary`: (`::Dict{Int64, Dict{String, Any}}`) event library

# Returns
- `events`: (`::Dict{Int64, Dict{String, Any}}`) event library
"""
function read_events(io, scale; type=nothing, format="%i "*"%f "^(length(scale)), eventLibrary=Dict())
    eventLength = length(scale) + 1
    fmt = Scanf.Format(format)
    args = Tuple([f == "%i" ? Int : (f == "%f" ? Float64 : (f == "%c" ? Char : String)) for f in split(format)])
    while true
        line = readline(io)
        isempty(line) && break
        r, data... = scanf(line, fmt, args...)
        r == eventLength || break #Break if not all values read
        id = floor(Int, data[1])
        data = [d isa Number ? s*d : d for (s, d) in zip(scale, data[2:end])]
        eventLibrary[id] = type === nothing ? Dict("data"=>data) : Dict("data"=>data, "type"=>type)
    end
    return eventLibrary
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
    extensionSpecLibrary[ext_id] = read_events(io, KomaMRIBase.get_scale(ext_type); format="%i "*KomaMRIBase.get_scanf_format(ext_type))
end
function read_extensions(io, ext_string, ext_type::Nothing, ext_id, extensionTypeLibrary, extensionSpecLibrary, required_extensions)
    if ext_string in required_extensions
        @error "Extension $ext_string is required by the sequence but not supported by KomaMRI reader"
    else
        @warn "Ignoring unsupported extension: $ext_string"
    end
end

"""
read_shapes Read the [SHAPES] section of a sequence file.
   library=read_shapes(fid) Read shapes from file identifier of an
   open MR sequence file and return a library of shapes.
"""
function read_shapes(io, forceConvertUncompressed)
    shapeLibrary = Dict()
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
    fix_first_last_grads!(obj::Dict)

Updates the `obj` dictionary with new first and last points for gradients.

# Notes:
- This function is "replicating" the following MATLAB code:
https://github.com/pulseq/pulseq/blob/v1.5.1/matlab/%2Bmr/%40Sequence/read.m#L325-L413
- We are updating the `gradLibrary` entries with the new first and last points, making them compatible with the v1.5.x format.
"""
function fix_first_last_grads!(obj::Dict, pulseq_version; simplify_shapes=true) 
    # Add first and last Pulseq points
    grad_prev_last = [0.0; 0.0; 0.0]
    for iB in 1:length(obj["blockEvents"])
        eventIDs = obj["blockEvents"][iB];
        block = get_block(obj, iB, pulseq_version; simplify_shapes=simplify_shapes)
        processedGradIDs = zeros(1, 3);    
        for iG in 1:3
            g_id = eventIDs[2+iG]
            g_id > 0 || continue
            g = obj["gradLibrary"][g_id]
            grad = g["data"]
            if g["type"] === 'g' # Arbitrary gradient waveform
                v1_5 = length(grad) == 6 # (1)amplitude (2)first_grads (3)last_grads (4)amp_shape_id (5)time_shape_id (6)delay 
                v1_4 = length(grad) == 4 # (1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
                v1_3 = length(grad) == 3 # (1)amplitude (2)amp_shape_id (3)delay
                v1_5 && continue # Already-updated (to v1.5 format) gradLibrary entry
                # From here, we are dealing with version == 1.5.x
                grad = v1_4 ? [grad[1]; grad_prev_last[iG]; 0.0; grad[2:end]] : [grad[1]; grad_prev_last[iG]; 0.0; grad[2]; 0; grad[3]]
                waveform = grad[1] * decompress_shape(obj["shapeLibrary"][grad[4]]...)
                if grad[6] > 0 # delay > 0
                    grad_prev_last[iG] = 0.0 
                end
                if grad[5] != 0 # time-shaped case
                    grad[3] = waveform[end]
                else # uniformly-shaped case
                    odd_step1 = [grad[2]; 2 * waveform]
                    odd_step2 = odd_step1 .* (mod.(1:length(odd_step1), 2) * 2 .- 1)
                    waveform_odd_rest = cumsum(odd_step2) .* (mod.(1:length(odd_step2), 2) * 2 .- 1)
                    grad[3] = waveform_odd_rest[end]
                end
                grad_prev_last[iG] = dur(block.GR[iG]) + eps(Float64) < dur(block) ? 0 : grad[3] # Bookkeeping for the next gradient
                if iG>1 && any(processedGradIDs[1:iG] .== g_id)
                    continue # avoid repeated updates if the same gradient is applied on differen gradient axes
                end
                processedGradIDs[iG] = g_id;
                obj["gradLibrary"][g_id]["data"] = grad # Update the gradLibrary entry
            else # Trapezoidal gradient waveform
                grad_prev_last[iG] = 0.0
            end
        end
    end
end

"""
    amp_shape, time_shape = simplify_waveforms(amp_shape, time_shape)

Simplifies the amplitude and time waveforms to a single value if they are constant.
"""
function simplify_waveforms(amp_shape, time_shape)
    if all(x->x==amp_shape[1], amp_shape)
        amp_shape = amp_shape[1]
        time_shape = sum(time_shape)
    elseif all(x->x==time_shape[1], time_shape)
        amp_shape = amp_shape
        time_shape = sum(time_shape)
    end
    return amp_shape, time_shape
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
function read_seq(filename; simplify_shapes=true)
    @info "Loading sequence $(basename(filename)) ..."
    pulseq_version = v"0.0.0"
    gradLibrary = Dict()
    def = Dict()
    signature = nothing
    blockEvents = Dict()
    rfLibrary = Dict()
    adcLibrary = Dict()
    tmp_delayLibrary = Dict()
    shapeLibrary = Dict()
    extensionInstanceLibrary = Dict()
    extensionTypeLibrary = Dict()
    extensionSpecLibrary = Dict()
    #Reading file and storing data
    open(filename) do io
        while !eof(io)
            section = readline(io)
            if typeof(section) == String && (isempty(section) || section[1] == '#')
                #skip useless line
            elseif section == "[DEFINITIONS]"
                def = read_definitions(io)
            elseif  section == "[VERSION]"
                pulseq_version = read_version(io)
            elseif  section == "[BLOCKS]"
                if pulseq_version == v"0.0.0"
                    @error "Pulseq file MUST include [VERSION] section prior to [BLOCKS] section"
                end
                blockEvents = read_blocks(io, pulseq_version)
            elseif  section == "[RF]"
                if pulseq_version >= v"1.5.0"
                    rfLibrary = read_events(io, [1/γ 1 1 1 1e-6 1e-6 1 1 1 1 1]; format="%i "*"%f "^(10)*"%c ") # this is 1.5.x format
                elseif pulseq_version >= v"1.4.0"
                    rfLibrary = read_events(io, [1/γ 1 1 1 1e-6 1 1]) # this is 1.4.x format
                else
                    rfLibrary = read_events(io, [1/γ 1 1 1e-6 1 1]) # this is 1.3.x and below
                    # we will have to scan through the library later after all the shapes have been loaded
                end
            elseif  section == "[GRADIENTS]"
                if pulseq_version >= v"1.5.0"
                    gradLibrary = read_events(io, [1/γ 1/γ 1/γ 1 1 1e-6]; type='g', eventLibrary=gradLibrary) # this is 1.5.x format
                elseif pulseq_version >= v"1.4.0"
                    gradLibrary = read_events(io, [1/γ 1 1 1e-6]; type='g', eventLibrary=gradLibrary) # this is 1.4.x format
                else
                    gradLibrary = read_events(io, [1/γ 1 1e-6]; type='g', eventLibrary=gradLibrary) # this is 1.3.x and below
                end
            elseif  section == "[TRAP]"
                gradLibrary = read_events(io, [1/γ 1e-6 1e-6 1e-6 1e-6]; type='t', eventLibrary=gradLibrary);
            elseif  section == "[ADC]"
                if pulseq_version >= v"1.5.0"
                    adcLibrary = read_events(io, [1 1e-9 1e-6 1 1 1 1 1]) # this is 1.5.x format
                else
                    adcLibrary = read_events(io, [1 1e-9 1e-6 1 1]) # this is 1.4.x and below
                end
            elseif  section == "[DELAYS]"
                if pulseq_version >= v"1.4.0"
                    @error "Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section"
                end
                tmp_delayLibrary = read_events(io, 1e-6);
            elseif  section == "[SHAPES]"
                shapeLibrary = read_shapes(io, (pulseq_version.major == 1 && pulseq_version.minor < 4))
            elseif  section == "[EXTENSIONS]"
                extensionInstanceLibrary = read_events(io, [1 1 1]; format="%i "*"%i "^(3), eventLibrary=extensionInstanceLibrary)
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
        # scan through the gradient objects and update 't'-s (trapezoids) und 'g'-s (free-shape gradients)
        for i in keys(gradLibrary)
            if gradLibrary[i]["type"] == 't'
                #(1)amplitude (2)rise (2)flat (3)fall (4)delay
                if gradLibrary[i]["data"][2] == 0 #rise
                    if abs(gradLibrary[i]["data"][1]) == 0 && gradLibrary[i]["data"][3] > 0
                        gradLibrary[i]["data"][3] -= def["gradRasterTime"]
                        gradLibrary[i]["data"][2]  = def["gradRasterTime"]
                    end
                end
                if gradLibrary[i]["data"][4] == 0 #delay
                    if abs(gradLibrary[i]["data"][1]) == 0 && gradLibrary[i]["data"][3] > 0
                        gradLibrary[i]["data"][3] -= def["gradRasterTime"]
                        gradLibrary[i]["data"][4]  = def["gradRasterTime"]
                    end
                end
            end
        end
    end
    verify_signature!(filename, signature; pulseq_version=pulseq_version)
    isempty(def) && (def = DEFAULT_DEFINITIONS)
    #Sequence
    obj = Dict(
        "blockEvents"=>blockEvents,
        "gradLibrary"=>gradLibrary,
        "rfLibrary"=>rfLibrary,
        "adcLibrary"=>adcLibrary,
        "tmp_delayLibrary"=>tmp_delayLibrary,
        "shapeLibrary"=>shapeLibrary,
        "extensionInstanceLibrary"=>extensionInstanceLibrary,
        "extensionTypeLibrary"=>extensionTypeLibrary,
        "extensionSpecLibrary"=>extensionSpecLibrary,
        "definitions"=>def
    )
    # Add first and last points for gradients #320 for version <= 1.4.2
    if pulseq_version < v"1.5.0"
        fix_first_last_grads!(obj, pulseq_version; simplify_shapes=simplify_shapes)
    end
    #Transforming Dictionary to Sequence object
    #This should only work for Pulseq files >=1.4.0
    seq = Sequence()
    for i = 1:length(blockEvents)
        seq += get_block(obj, i, pulseq_version; simplify_shapes=simplify_shapes)
    end
    # Final details
    #Temporary hack
    seq.DEF = merge(obj["definitions"], seq.DEF)
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
- `gradLibrary`: (`::Dict{K, V}`) the "gradLibrary" dictionary
- `shapeLibrary`: (`::Dict{K, V}`) the "shapeLibrary" dictionary
- `Δt_gr`: (`::Float64`, `[s]`) gradient raster time
- `i`: (`::Int64`) index of the axis in the block event

# Returns
- `grad`: (::Grad) Gradient struct
"""
function get_Grad(gradLibrary, shapeLibrary, Δt_gr, i; simplify_shapes=true)
    gr = Grad(0.0,0.0)
    haskey(gradLibrary, i) || return gr
    if gradLibrary[i]["type"] === 't' # Trapezoidal gradient waveform
        #(1)amplitude (2)rise (3)flat (4)fall (5)delay
        g_A, g_rise, g_T, g_fall, g_delay = gradLibrary[i]["data"]
        gr = Grad(g_A,g_T,g_rise,g_fall,g_delay,0.0,0.0)
    elseif gradLibrary[i]["type"] === 'g' # Arbitrary gradient waveform
        g = gradLibrary[i]["data"]
        v1_5 = length(g) == 6 # (1)amplitude (2)first_grads (3)last_grads (4)amp_shape_id (5)time_shape_id (6)delay 
        v1_4 = length(g) == 4 # (1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
        v1_3 = length(g) == 3 # (1)amplitude (2)amp_shape_id (3)delay
        amplitude     =  g[1]
        first_grads   =  v1_5 ? g[2] : 0.0
        last_grads    =  v1_5 ? g[3] : 0.0
        amp_shape_id  = (v1_5 ? g[4] : g[2]) |> x->floor(Int64,x)
        time_shape_id = (v1_5 ? g[5] : (v1_4 ? g[3] : 0)) |> x->floor(Int64,x)
        delay         =  v1_5 ? g[6] : (v1_4 ? g[4] : g[3])
        #Amplitude
        gA = amplitude * decompress_shape(shapeLibrary[amp_shape_id]...)
        Ngr = length(gA) - 1
        #Creating timings
        if time_shape_id == 0 #no time waveform. Default time raster
            gT = Ngr * Δt_gr
            rise, fall = Δt_gr/2, Δt_gr/2
        elseif time_shape_id == -1 #New in pulseq 1.5.x: no time waveform. 1/2 of the default time raster
            gT = Ngr * Δt_gr / 2
            rise, fall = Δt_gr/2, Δt_gr/2
        else
            gt = decompress_shape(shapeLibrary[time_shape_id]...)
            gt[1] == 0 || @warn "Gradient time shape $time_shape_id starting at a non-zero value $(gt[1]). This is not recommended and may not be supported properly\n (see https://github.com/pulseq/pulseq/issues/188#issuecomment-3541588756) " maxlog=1
            delay += gt[1] * Δt_gr # offset due to the shape starting at a non-zero value. This case 
            gT = diff(gt) * Δt_gr
            rise, fall = 0.0, 0.0
        end
        if simplify_shapes gA, gT = simplify_waveforms(gA, gT) end
        gr = Grad(gA, gT, rise, fall, delay, first_grads, last_grads)
    end
    return gr
end

"""
    rf = get_RF(rfLibrary, shapeLibrary, Δt_rf, i)

Reads the RF. It is used internally by [`get_block`](@ref).

# Arguments
- `rfLibrary`: (`::Dict{K, V}`) the "rfLibrary" dictionary
- `shapeLibrary`: (`::Dict{K, V}`) the "shapeLibrary" dictionary
- `Δt_rf`: (`::Float64`, `[s]`) RF raster time
- `i`: (`::Int64`) index of the RF in the block event

# Returns
- `rf`: (`1x1 ::Matrix{RF}`) RF struct
"""
function get_RF(rfLibrary, shapeLibrary, Δt_rf, i; simplify_shapes=true)
    rf = RF(0.0,0.0)
    haskey(rfLibrary, i) || return [rf;;], false
    #Unpacking
    r = rfLibrary[i]["data"]
    v1_5 = length(r) == 11 # (1)amplitude (2)mag_id (3)phase_id (4)time_shape_id (5)center (6)delay (7)freq_ppm (8)phase_ppm (9)freq (10)phase (11)use
    v1_4 = length(r) == 7  # (1)amplitude (2)mag_id (3)phase_id (4)time_shape_id (5)delay (6)freq (7)phase
    v1_3 = length(r) == 6  # (1)amplitude (2)mag_id (3)phase_id (4)delay (5)freq (6)phase
    amplitude     =  r[1]
    mag_id        =  r[2] |> x->floor(Int64,x)
    phase_id      =  r[3] |> x->floor(Int64,x)
    time_shape_id =  v1_3 ? 0 : (r[4] |> x->floor(Int64,x)); add_half_Δt_rf = !(time_shape_id>0)
    center        =  v1_5 ? r[5]  : 0.0
    delay         = (v1_5 ? r[6]  : (v1_4 ? r[5] : r[4])) + (add_half_Δt_rf) * Δt_rf/2
    freq_ppm      =  v1_5 ? r[7]  : 0.0
    phase_ppm     =  v1_5 ? r[8]  : 0.0
    freq          =  v1_5 ? r[9]  : (v1_4 ? r[6] : r[5])
    phase         =  v1_5 ? r[10] : (v1_4 ? r[7] : r[6])
    use           =  v1_5 ? r[11] : 'u'
    #Amplitude and phase waveforms
    if amplitude != 0 && mag_id != 0
        rfA = decompress_shape(shapeLibrary[mag_id]...)
        rfϕ = decompress_shape(shapeLibrary[phase_id]...)
        @assert all(rfϕ.>=0) "[RF id $i] Phase waveform rfϕ must have non-negative samples (1.>=rfϕ.>=0). "
        Nrf = shapeLibrary[mag_id][1] - 1
        rfAϕ = amplitude .* rfA .* exp.(1im*(2π*rfϕ .+ phase))
    else
        rfAϕ = ComplexF64[0.0]
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
    if simplify_shapes rfAϕ, rfT = simplify_waveforms(rfAϕ, rfT) end
    use = KomaMRIBase.get_RF_use_from_char(Val(use))
    if v1_5 # for version 1.5.x
        rf = RF(rfAϕ,rfT,freq,delay,center,use)
    else # for version 1.4.x and below
        rf = RF(rfAϕ,rfT,freq,delay)
    end
    return [rf;;], add_half_Δt_rf
end

"""
    adc = get_ADC(adcLibrary, i)

Reads the ADC. It is used internally by [`get_block`](@ref).

# Arguments
- `adcLibrary`: (`::Dict{String, Any}`) the "adcLibrary" dictionary
- `i`: (`::Int64`) index of the adc in the block event

# Returns
- `adc`: (`1x1 ::Vector{ADC}`) ADC struct
"""
function get_ADC(adcLibrary, i)
    adc = ADC(0, 0)
    haskey(adcLibrary, i) || return [adc], 0.0
    #Unpacking
    a = adcLibrary[i]["data"]
    v1_5 = length(a) == 8 # (1)num (2)dwell (3)delay (4)freq_ppm (5)phase_ppm (6)freq (7)phase (8)phase_shape_id
    v1_4 = length(a) == 5 # (1)num (2)dwell (3)delay (4)freq (5)phase
    num       = a[1] |> x->floor(Int64,x)
    dwell     = a[2]
    delay     = a[3] + dwell/2
    freq_ppm  = v1_5 ? a[4] : 0.0
    phase_ppm = v1_5 ? a[5] : 0.0
    freq      = v1_5 ? a[6] : a[4]
    phase     = v1_5 ? a[7] : a[5]
    phase_id  = v1_5 ? a[8] : 0
    #Definition
    T = (num-1) * dwell
    adc = ADC(num,T,delay,freq,phase)
    return [adc], delay + T - dwell/2
end

"""
    seq = get_block(obj, i, pulseq_version; simplify_shapes=true)

Block sequence definition. Used internally by [`read_seq`](@ref).

# Arguments
- `obj`: (`::Dict{String, Any}`) main dictionary
- `i`: (`::Int64`) index of a block event
- `pulseq_version`: (`::VersionNumber`) Pulseq version

# Keywords
- `simplify_shapes`: (`::Bool`) whether to simplify the shapes of the gradients and RFs

# Returns
- `s`: (`::Sequence`) block Sequence struct
"""
function get_block(obj, i, pulseq_version; simplify_shapes=true)
    #Unpacking
    idur, irf, igx, igy, igz, iadc, iext = obj["blockEvents"][i]
    #Gradient definition
    Δt_gr = obj["definitions"]["GradientRasterTime"]
    Gx = get_Grad(obj["gradLibrary"], obj["shapeLibrary"], Δt_gr, igx; simplify_shapes=simplify_shapes)
    Gy = get_Grad(obj["gradLibrary"], obj["shapeLibrary"], Δt_gr, igy; simplify_shapes=simplify_shapes)
    Gz = get_Grad(obj["gradLibrary"], obj["shapeLibrary"], Δt_gr, igz; simplify_shapes=simplify_shapes)
    G = reshape([Gx;Gy;Gz],3,1) #[Gx;Gy;Gz;;]
    #RF definition
    Δt_rf = obj["definitions"]["RadiofrequencyRasterTime"]
    R, add_half_Δt_rf = get_RF(obj["rfLibrary"], obj["shapeLibrary"], Δt_rf, irf; simplify_shapes=simplify_shapes)
    #ADC definition
    A, adc_dur = get_ADC(obj["adcLibrary"], iadc)
    #DUR
    max_dur = max(dur(Gx), dur(Gy), dur(Gz), dur(R[1]) + (add_half_Δt_rf) * Δt_rf/2, adc_dur)
    if pulseq_version >= v"1.4.0" # Explicit block duration (in units of blockDurationRaster)
        duration = idur * obj["definitions"]["BlockDurationRaster"]
        @assert duration ≈ max_dur || duration >= max_dur "Block duration must be greater than or approximately equal to the duration of the block events"
        D = Float64[duration]
    else # Block duration as the maximum between the delay and the duration of the block events
        idelay = idur > 0 ? obj["tmp_delayLibrary"][idur]["data"][1] : 0.0
        D = Float64[max(idelay, max_dur)]
    end
    #Extensions
    E = get_extension(obj["extensionInstanceLibrary"], obj["extensionTypeLibrary"], obj["extensionSpecLibrary"], iext)
    # Definitition
    DEF = Dict{String,Any}()
    #Sequence block definition
    return Sequence(G,R,A,D,E,DEF)
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
    EXT = [Extension[]]
    haskey(extensionInstanceLibrary, i) || return EXT
    # type ref next_id
    # next_id of 0 terminates the list
    type, ref, next_id = extensionInstanceLibrary[i]["data"]
    while true
        length(extensionTypeLibrary) < type ? (@warn "extension type n°$type does not exist"; break) : nothing
        push!(EXT[1], extensionTypeLibrary[type](extensionSpecLibrary[type][ref]["data"]...))
        if next_id == 0
            break
        else
            type, ref, next_id = extensionInstanceLibrary[next_id]["data"]
        end
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
Base.@kwdef struct PulseqExportContext
    seq::Sequence
    filename::String
    version::VersionNumber = v"1.5.1"
    blockDurationRaster::Float64 = get_blockDurationRaster(seq)
    gradientRasterTime::Float64 = get_gradientRasterTime(seq)
    rfRasterTime::Float64 = get_rfRasterTime(seq)
    adcRasterTime::Float64 = get_adcRasterTime(seq)
end

get_blockDurationRaster(seq::Sequence) = get(seq.DEF, "BlockDurationRaster", DEFAULT_DEFINITIONS["BlockDurationRaster"])
get_gradientRasterTime(seq::Sequence)  = get(seq.DEF, "GradientRasterTime", DEFAULT_DEFINITIONS["GradientRasterTime"])
get_rfRasterTime(seq::Sequence)        = get(seq.DEF, "RadiofrequencyRasterTime", DEFAULT_DEFINITIONS["RadiofrequencyRasterTime"])
get_adcRasterTime(seq::Sequence)       = get(seq.DEF, "AdcRasterTime", DEFAULT_DEFINITIONS["AdcRasterTime"])

Base.@kwdef struct PulseqBlock
    id::Int
    duration::Int = 0
    rf::Int = 0
    gx::Int = 0
    gy::Int = 0
    gz::Int = 0
    adc::Int = 0
    ext::Int = 0
end

# - Dict keys are the ids of the objects
# - Dict values are tuples containing the rest of the values of the blocks/events/shapes
Base.@kwdef struct PulseqExportAssets
    blocks::Dict{Int,PulseqBlock} = Dict{Int,PulseqBlock}()
    rf::Dict{Int,Tuple{Float64,Int,Int,Int,Float64,Int,Float64,Float64,Float64,Float64,Char}} = Dict{Int,Tuple{Float64,Int,Int,Int,Float64,Int,Float64,Float64,Float64,Float64,Char}}()
    gradients::Dict{Int,Tuple{Float64,Float64,Float64,Int,Int,Int}} = Dict{Int,Tuple{Float64,Float64,Float64,Int,Int,Int}}()
    trapezoids::Dict{Int,Tuple{Float64,Int,Int,Int,Int}} = Dict{Int,Tuple{Float64,Int,Int,Int,Int}}()
    adc::Dict{Int,Tuple{Int,Float64,Int,Float64,Float64,Float64,Float64,Int}} = Dict{Int,Tuple{Int,Float64,Int,Float64,Float64,Float64,Float64,Int}}()
    shapes::Dict{Int,Tuple{Int,Vector{Float64}}} = Dict{Int,Tuple{Int,Vector{Float64}}}() # shapes are stored compressed if selected
    extensionInstances::Dict{Int,Tuple{Int,Int,Int}} = Dict{Int,Tuple{Int,Int,Int}}()
    extensionTypes::Dict{Int,Type{<:Extension}} = Dict{Int,Type{<:Extension}}()
    extensionSpecs::Dict{Int,Dict{Int,Extension}} = Dict{Int,Dict{Int,Extension}}()
end

"""
    collect_pulseq_assets(ctx::PulseqExportContext) -> PulseqExportAssets

Create the Pulseq export dictionaries required to serialize `ctx.seq` into the 1.5.1 file format.
This function is responsible for deduplicating reusable objects (RF, gradients, shapes, etc.)
and for translating each sequence block into the integer lookups (i.e., the assets) expected by the specification.
"""
function collect_pulseq_assets(ctx::PulseqExportContext)
    assets = PulseqExportAssets()
    seq = ctx.seq
    
    # First step: Collect all unique extension vectors and register them once
    # This ensures that blocks sharing the same extensions reuse the same instance IDs
    # We use parallel arrays to track extension vectors and their IDs (since vectors can't be dict keys)
    extension_vectors = Vector{Extension}[]
    extension_vector_ids = Int[]
    
    for (idx, block) in enumerate(seq)
        ext_vec = block.EXT[1]
        if length(ext_vec) > 0
            # Check if this extension vector has already been registered by comparing content
            matching_idx = findfirst(registered_ext -> 
                length(registered_ext) == length(ext_vec) && 
                all(e1 ≈ e2 for (e1, e2) in zip(registered_ext, ext_vec)),
                extension_vectors
            )
            if matching_idx === nothing
                # Register this extension vector
                ext_id = register_ext!(assets, ext_vec, ctx)
                push!(extension_vectors, ext_vec)
                push!(extension_vector_ids, ext_id)
            end
        end
    end
    
    # Second step: Process all blocks and use pre-registered extension IDs
    for (idx, block) in enumerate(seq)
        # 1. RF: deduplicate waveform in `assets.rf` and store the ID
        rf = block.RF[1]
        rf_id = register_rf!(assets, rf.A, rf.T, rf.Δf, rf.delay, rf.center,rf.use, ctx)

        # 2. Gradients: do the same per axis, store ids `gx`, `gy`, `gz`
        gx, gy, gz = block.GR
        gx_id = register_grad!(assets, gx.A, gx.T, gx.rise, gx.fall, gx.delay, gx.first, gx.last, ctx)
        gy_id = register_grad!(assets, gy.A, gy.T, gy.rise, gy.fall, gy.delay, gy.first, gy.last, ctx)
        gz_id = register_grad!(assets, gz.A, gz.T, gz.rise, gz.fall, gz.delay, gz.first, gz.last, ctx)

        # 3. ADC: deduplicate and store the ID
        adc = block.ADC[1]
        adc_id = register_adc!(assets, adc.N, adc.T, adc.delay, adc.Δf, adc.ϕ, ctx)

        # 4. Calculate block duration: must be >= max duration of all events
        # After quantization, event durations may have changed, so we need to recalculate
        max_event_duration = max(dur(gx), dur(gy), dur(gz), dur(rf), dur(adc))
        block_duration = max(dur(block), max_event_duration)
        # Round to block duration raster (round up to ensure >= max_event_duration)
        duration = ceil(Int, block_duration / ctx.blockDurationRaster)

        # 5. Extensions: find pre-registered ID by content comparison
        ext_vec = block.EXT[1]
        if length(ext_vec) > 0
            matching_idx = findfirst(registered_ext -> 
                length(registered_ext) == length(ext_vec) && 
                all(e1 ≈ e2 for (e1, e2) in zip(registered_ext, ext_vec)),
                extension_vectors
            )
            ext_id = matching_idx === nothing ? 0 : extension_vector_ids[matching_idx]
        else
            ext_id = 0
        end

        assets.blocks[idx] = PulseqBlock(idx, duration, rf_id, gx_id, gy_id, gz_id, adc_id, ext_id)
    end
    return assets
end

"""
    id = register_rf!(assets, A, T, Δf, delay, use, ctx)
"""
function register_rf!(assets::PulseqExportAssets, A, T, Δf, delay, center, use, ctx::PulseqExportContext)
    iszero(maximum(abs.(A))) && return 0
    Δt_rf = ctx.rfRasterTime
    mag_id, phase_id, time_id = register_rf_shapes!(assets.shapes, A, T, Δt_rf; compress=true)
    amp          = γ * maximum(abs.(A)) # from T to Hz (nucleus-dependent)
    center       = center * 1e6 # from s to us         # Would get_RF_center(RF(A, T, Δf, delay, center, use)) be enough?
    delay        = round(Int, (delay - (time_id==0) * Δt_rf/2) * 1e6) # from s to us
    freq_ppm     = 0.0
    phase_ppm    = 0.0
    freq_offset  = Δf
    phase_offset = 0.0
    use          = KomaMRIBase.get_char_from_RF_use(use)
    aux = (amp, mag_id, phase_id, time_id, center, delay, freq_ppm, phase_ppm, freq_offset, phase_offset, use)
    return _store_event!(assets.rf, aux)
end

# A and T are numbers (pulse waveform)
function register_rf_shapes!(shapes, A, T, Δrf; compress = true)
    mag_id    = _store_shape!(shapes, [1.0, 1.0]; compress=compress)
    phase_id  = _store_shape!(shapes, [0.0, 0.0]; compress=compress)
    time_id   = _store_shape!(shapes, [0.0, T] ./ Δrf; compress=compress)
    return mag_id, phase_id, time_id
end
# A is a vector (uniformly-sampled waveform)
function register_rf_shapes!(shapes, A::Vector, T, Δrf; compress = true)
    n_samples = length(A)
    mag_id    = _store_shape!(shapes, abs.(A) ./ maximum(abs.(A)); compress=compress)
    phase_id  = _store_shape!(shapes, mod.(angle.(A) / (2π), 1.0); compress=compress)
    t_vector  = collect(range(0, T, length=n_samples)) ./ Δrf
    time_id   = _store_shape!(shapes, t_vector; compress=compress)
    return mag_id, phase_id, time_id
end
# A and T are vectors (time-shaped waveform)
function register_rf_shapes!(shapes, A::Vector, T::Vector, Δrf; compress = true)
    mag_id    = _store_shape!(shapes, abs.(A) ./ maximum(abs.(A)); compress=compress)
    phase_id  = _store_shape!(shapes, mod.(angle.(A) / (2π), 1.0); compress=compress)
    t_vector = cumsum(vcat(0.0, T ./ Δrf))
    time_id   = _store_shape!(shapes, t_vector; compress=compress)
    return mag_id, phase_id, time_id
end

"""
    id = register_grad!(assets, A, T, rise, fall, delay, first, last, ctx)
"""
# Arbitrary gradient waveform (into [GRADIENTS] section) (rise and fall ARE NOT USED)
function register_grad!(assets::PulseqExportAssets, A::Vector, T, rise, fall, delay, first, last, ctx::PulseqExportContext)
    iszero(maximum(abs.(A))) && return 0
    shape_id, time_id = register_grad_shapes!(assets.shapes, A, T, ctx.gradientRasterTime; compress=true)
    amp   = γ * maximum(abs.(A)) # from T/m to Hz/m (nucleus-dependent)
    delay = round(Int, delay * 1e6) # from s to us
    first = γ * first # from T/m to Hz/m 
    last  = γ * last  # from T/m to Hz/m
    aux = (amp, first, last, shape_id, time_id, delay)
    return _store_event!(assets.gradients, aux)
end
# Trapezoidal gradient waveform (into [TRAP] section) (first and last ARE NOT USED)
function register_grad!(assets::PulseqExportAssets, A::Real, T, rise, fall, delay, first, last, ctx::PulseqExportContext)
    iszero(A) && return 0
    amp   = γ * A # from T/m to Hz/m (nucleus-dependent)  
    rise  = round(Int, rise  * 1e6) # from s to us
    flat  = round(Int, T     * 1e6) # from s to us
    fall  = round(Int, fall  * 1e6) # from s to us
    delay = round(Int, delay * 1e6) # from s to us
    aux = (amp, rise, flat, fall, delay)
    return _store_event!(assets.trapezoids, aux)
end 

# A is a vector and T is a number (uniformly-sampled waveform)
function register_grad_shapes!(shapes, A::Vector, Δgr; compress = true)
    shape_id  = _store_shape!(shapes, A ./ maximum(abs.(A)); compress=compress)
    return shape_id, 0
end
# A and T are vectors (time-shaped waveform)
function register_grad_shapes!(shapes, A::Vector, T::Vector, Δgr; compress = true)
    shape_id = _store_shape!(shapes, A ./ maximum(abs.(A)); compress=compress)
    t_vector = cumsum([0;  T]) ./ Δgr
    time_id  = _store_shape!(shapes, t_vector; compress=compress)
    return shape_id, time_id
end

"""
    id = register_adc!(assets, N, T, delay, Δf, ϕ, ctx)
"""
function register_adc!(assets::PulseqExportAssets, N, T, delay, Δf, ϕ, ctx::PulseqExportContext)
    iszero(N) && return 0
    dwell_s = N == 1 ? T : T / (N - 1)
    dwell = dwell_s * 1e9 # from s to ns
    del = round(Int, (delay - dwell_s/2) * 1e6) # from s to us, subtract dwell/2 
    freq_ppm = 0.0
    phase_ppm = 0.0
    freq = Δf
    phase = ϕ
    phase_id = 0 # TODO: implement phase shape in Koma
    aux = (N, dwell, del, freq_ppm, phase_ppm, freq, phase, phase_id)
    return _store_event!(assets.adc, aux)
end

"""
    id = register_ext!(assets, ext, ctx)
"""
function register_ext!(assets::PulseqExportAssets, ext::Vector{Extension}, ctx::PulseqExportContext)
    length(ext) == 0 && return 0
    instance_ids = Int[]
    for e in ext
        ext_id      = _store_event!(assets.extensionTypes, typeof(e))
        
        if !haskey(assets.extensionSpecs, ext_id)
            assets.extensionSpecs[ext_id] = Dict{Int,Extension}()
        end
        
        ref         = _store_event!(assets.extensionSpecs[ext_id], e)
        instance_id = _store_event!(assets.extensionInstances, (ext_id, ref, 0))
        push!(instance_ids, instance_id)
    end
    
    for i in 1:(length(instance_ids) - 1)
        current_id = instance_ids[i]
        next_id = instance_ids[i + 1]
        
        # Get the current tuple and update it with the correct next_id
        type_id, ref_id, _ = assets.extensionInstances[current_id]
        assets.extensionInstances[current_id] = (type_id, ref_id, next_id)
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
)
    payload = compress ? compress_shape(samples) : (length(samples), samples)
    for (k, existing) in shapes
        if existing == payload
            return k
        end
    end
    new_id = maximum(keys(shapes); init=0) + 1
    shapes[new_id] = payload
    return new_id
end

"""
    emit_pulseq(io::IO, ctx::PulseqExportContext, assets::PulseqExportAssets)

Write the Pulseq sections to `io`, using the already prepared `assets` and `ctx` 
"""
function emit_pulseq(io::IO, ctx::PulseqExportContext, assets::PulseqExportAssets)
    emit_header_comment!(io)
    emit_version_section!(io, ctx)
    emit_definitions_section!(io, ctx)
    if isempty(assets.blocks)
        @warn "Pulseq export: block table is empty; populate it via `collect_pulseq_assets`."
    else
        emit_blocks_section!(io, assets)
    end
    if !isempty(assets.rf)
        emit_rf_section!(io, assets)
    end
    if !isempty(assets.gradients)
        emit_gradients_section!(io, assets)
    end
    if !isempty(assets.trapezoids)
        emit_trap_section!(io, assets)
    end
    if !isempty(assets.adc)
        emit_adc_section!(io, assets)
    end
    if !isempty(assets.extensionInstances)
        emit_extension_section!(io, assets)
    end
    if !isempty(assets.shapes)
        emit_shapes_section!(io, assets)
    end
end

emit_header_comment!(io::IO) = write(io, "# Pulseq sequence file\n# Created by KomaMRI.jl\n\n")

function emit_version_section!(io::IO, ctx::PulseqExportContext)
    write(io, "[VERSION]\n")
    write(io, "major $(ctx.version.major)\n")
    write(io, "minor $(ctx.version.minor)\n")
    write(io, "revision $(ctx.version.patch)\n\n")
end

function emit_definitions_section!(io::IO, ctx::PulseqExportContext)
    write(io, "[DEFINITIONS]\n")
    write(io, "BlockDurationRaster $(ctx.blockDurationRaster)\n")
    write(io, "GradientRasterTime $(ctx.gradientRasterTime)\n")
    write(io, "RadiofrequencyRasterTime $(ctx.rfRasterTime)\n")
    write(io, "AdcRasterTime $(ctx.adcRasterTime)\n")
    for (key, value) in ctx.seq.DEF
        write(io, "$key $value\n")
    end
    write(io, "\n")
end

function emit_blocks_section!(io::IO, assets::PulseqExportAssets)
    write(io, "# Format of blocks:\n")
    write(io, "# NUM DUR RF  GX  GY  GZ  ADC  EXT\n")
    write(io, "[BLOCKS]\n")
    isempty(assets.blocks) && return
    block_ids = sort(collect(keys(assets.blocks)))
    for id in block_ids
        block = assets.blocks[id]
        values = [id, block.duration, block.rf, block.gx, block.gy, block.gz, block.adc, block.ext]
        max_lengths = zeros(Int, length(values))
        for (i, val) in enumerate(values)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, values, max_lengths)
    end
end

function emit_rf_section!(io::IO, assets::PulseqExportAssets)
    write(io, "\n# Format of RF events:\n")
    write(io, "# id amp mag_id phase_id time_id center delay freq_ppm phase_ppm freq_off phase_off use\n")
    write(io, "# ..  Hz     ..       ..       ..    us    us      ppm   rad/MHz       Hz       rad  ..\n")
    write(io, "# Field 'use' is the initial of:\n")
    write(io, "# excitation refocusing inversion saturation preparation other undefined\n")
    write(io, "[RF]\n")
    isempty(assets.rf) && return
    rf_ids = sort(collect(keys(assets.rf)))
    for id in rf_ids
        rf_data = assets.rf[id]
        values = [id, rf_data[1], rf_data[2], rf_data[3], rf_data[4], rf_data[5], rf_data[6], rf_data[7], rf_data[8], rf_data[9], rf_data[10], rf_data[11]]
        max_lengths = zeros(Int, 12)
        for (i, val) in enumerate(values)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, values, max_lengths)
    end
end

function emit_gradients_section!(io::IO, assets::PulseqExportAssets)
    write(io, "\n# Format of arbitrary gradient events:\n")
    write(io, "# id      amp      first      last  shape_id  time_id  delay\n")
    write(io, "# ..     Hz/m       Hz/m      Hz/m        ..       ..     us\n")
    write(io, "[GRADIENTS]\n")
    isempty(assets.gradients) && return
    grad_ids = sort(collect(keys(assets.gradients)))
    for id in grad_ids
        grad_data = assets.gradients[id]
        vals = (id, grad_data[1], grad_data[2], grad_data[3], grad_data[4], grad_data[5], grad_data[6])
        max_lengths = zeros(Int, 7)
        for (i, val) in enumerate(vals)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, vals, max_lengths)
    end
end

function emit_trap_section!(io::IO, assets::PulseqExportAssets)
    write(io, "\n# Format of trapezoid gradient events:\n")
    write(io, "# id      amp      rise  flat  fall  delay\n")
    write(io, "# ..     Hz/m        us    us    us     us\n")
    write(io, "[TRAP]\n")
    isempty(assets.trapezoids) && return
    trap_ids = sort(collect(keys(assets.trapezoids)))
    for id in trap_ids
        trap_data = assets.trapezoids[id]
        vals = (id, trap_data[1], trap_data[2], trap_data[3], trap_data[4], trap_data[5])
        max_lengths = zeros(Int, 6)
        for (i, val) in enumerate(vals)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, vals, max_lengths)
    end
end

function emit_adc_section!(io::IO, assets::PulseqExportAssets)
    write(io, "\n# Format of ADC events:\n")
    write(io, "# id  num  dwell  delay  freq_ppm  phase_ppm  freq  phase  phase_id\n")
    write(io, "# ..   ..     ns     us       ppm        ppm    Hz    rad        ..\n")
    write(io, "[ADC]\n")
    isempty(assets.adc) && return
    adc_ids = sort(collect(keys(assets.adc)))
    for id in adc_ids
        adc_data = assets.adc[id]
        vals = (id, adc_data[1], adc_data[2], adc_data[3], adc_data[4], adc_data[5], adc_data[6], adc_data[7], adc_data[8])
        max_lengths = zeros(Int, 9)
        for (i, val) in enumerate(vals)
            str = _format_value(val)
            max_lengths[i] = max(max_lengths[i], length(str))
        end
        _format_row(io, vals, max_lengths)
    end
end

function emit_extension_section!(io::IO, assets::PulseqExportAssets)
    write(io, "\n# Format of extension events:\n")
    write(io, "# id  type  ref  next_id\n")
    write(io, "# Extension list is followed by extension specifications\n")
    write(io, "[EXTENSIONS]\n")
    isempty(assets.extensionInstances) && return
    ext_ids = sort(collect(keys(assets.extensionInstances)))
    for id in ext_ids
        ext_data = assets.extensionInstances[id]
        vals = (id, ext_data[1], ext_data[2], ext_data[3])
        max_lengths = length.(digits.(vals))
        _format_row(io, vals, max_lengths)
    end
    write(io, "\n")
    type_ids = sort(collect(keys(assets.extensionTypes)))
    for id in type_ids
        write(io, KomaMRIBase.extension_type_header(assets.extensionTypes[id]))
        write(io, "extension $(string(KomaMRIBase.get_symbol_from_EXT_type(assets.extensionTypes[id]))) $id\n")
        spec_ids = sort(collect(keys(assets.extensionSpecs[id])))
        for spec_id in spec_ids
            v = assets.extensionSpecs[id][spec_id]
            field_values = Tuple(getfield(v, f) for f in fieldnames(typeof(v)))
            scaled_values = [d isa Char ? d : s*d for (s, d) in zip(KomaMRIBase.get_scale(typeof(v)), field_values)]     
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

function emit_shapes_section!(io::IO, assets::PulseqExportAssets)
    write(io, "\n# Sequence shapes\n")
    write(io, "[SHAPES]\n\n")
    for id in sort(collect(keys(assets.shapes)))
        num_samples, samples = assets.shapes[id]
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
    write_seq(seq, filename)
"""
function write_seq(
    seq::Sequence, filename::String;
    blockDurationRaster = get_blockDurationRaster(seq),
    gradientRasterTime = get_gradientRasterTime(seq),
    rfRasterTime = get_rfRasterTime(seq),
    adcRasterTime = get_adcRasterTime(seq),
    signatureAlgorithm::AbstractString = "md5"
)
    @info "Saving sequence to $(basename(filename)) ..."
    ctx = PulseqExportContext(seq, filename, v"1.5.1", blockDurationRaster, gradientRasterTime, rfRasterTime, adcRasterTime)
    assets = collect_pulseq_assets(ctx)
    buffer = IOBuffer()
    emit_pulseq(buffer, ctx, assets)
    payload = take!(buffer)
    signature_hash = supported_signature_digest(signatureAlgorithm, payload)
    open(filename, "w") do io
        write(io, payload)
        write(io, '\n')
        emit_signature_section!(io, signatureAlgorithm, signature_hash)
    end
    return nothing
end
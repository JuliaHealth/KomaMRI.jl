# # add_format(format"Pulseq", "# Pulseq sequence file", ".seq", [:Koma=>UUID("a9056882-1c2f-47b4-848d-dcb4a04f1994")])
# module Pulseq
# This file is a copy of the MATLAB Pulseq v1.4.0 read functions

"""
read_version Read the [VERSION] section of a sequence file.
   defs=read_version(fid) Read Pulseq version from file
   identifier of an open MR sequence file and return it
"""
function read_version(io)
    major =    @scanf(readline(io), "major %i", Int)[end]
    minor =    @scanf(readline(io), "minor %i", Int)[end]
    revision = @scanf(readline(io), "revision %i", Int)[end]

    version_combined = 1_000_000*major+1_000*minor+revision

    @assert major == 1 "Unsupported version_major $major"
    if     version_combined < 1002000
        @error "Unsupported version $major.$minor.$revision, only file format revision 1.2.0 and above are supported"
    elseif version_combined < 1003001
        @warn "Loading older Pulseq format file (version $major.$minor.$revision) some code may not function as expected"
    end

    major, minor, revision, version_combined
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
    if !haskey(def,"BlockDurationRaster")       def["BlockDurationRaster"] = 1e-5       end
    if !haskey(def,"GradientRasterTime")        def["GradientRasterTime"] = 1e-5        end
    if !haskey(def,"RadiofrequencyRasterTime")  def["RadiofrequencyRasterTime"] = 1e-6  end
    if !haskey(def,"AdcRasterTime")             def["AdcRasterTime"] = 1e-7             end
    def
end

"""
read_blocks Read the [BLOCKS] section of a sequence file.
   library=read_blocks(fid) Read blocks from file identifier of an
   open MR sequence file and return the event table.
"""
function read_blocks(io, blockDurationRaster, version_combined)
    eventTable = Dict()
    blockDurations = Dict()
    delayIDs_tmp = Dict()
    while true
        if version_combined <= 1002001
            NumberBlockEvents = 7
        else
            NumberBlockEvents = 8
        end

        fmt = Scanf.Format("%i "^NumberBlockEvents)
        r, blockEvents... = scanf(readline(io), fmt, zeros(Int,NumberBlockEvents)...)

        if blockEvents[1] != 0
            if version_combined <= 1002001
                eventTable[blockEvents[1]] = [0 blockEvents[3:end]... 0]
            else
                eventTable[blockEvents[1]] = [0 blockEvents[3:end]...]
            end

            if version_combined >= 1004000
                blockDurations[blockEvents[1]] = blockEvents[2]*blockDurationRaster
            else
                delayIDs_tmp[blockEvents[1]] = blockEvents[2]
            end
        end

        r == NumberBlockEvents || break #Break on white space  # fixme: pulseq doesn't technically forbid empty lines, does it? or comments... this would trip! although p8 or 9 of spec says "each subsequent line..."
    end
    sort(eventTable), sort(blockDurations), sort(delayIDs_tmp)  # this is not correct! see docstring comment!
end
"""
read_blocks Read the [BLOCKS] section of a sequence file.
   library=read_blocks(fid) Read blocks from file identifier of an
   open MR sequence file and return the event table.
   produces same output as read_blocks, but does not use scanf,
   rather reads everything as int
"""
function read_blocks_split_instead_of_scanf(io, blockDurationRaster, version_combined)
    eventTable = Dict()
    blockDurations = Dict()
    delayIDs_tmp = Dict()
    if version_combined <= 1002001
        NumberBlockEvents = 7
    else
        NumberBlockEvents = 8
    end
    while true
        blockEvents = [parse(Int, d) for d in eachsplit(readline(io))]
        if length(blockEvents) == 0
            # empty line
            break
        end
        if length(blockEvents) != NumberBlockEvents
            error("Expected $NumberBlockEvents events but got $(length(blockEvents)).")
        end
        if blockEvents[1] != 0
            if version_combined <= 1002001
                eventTable[blockEvents[1]] = [0 blockEvents[3:end]... 0]
            else
                eventTable[blockEvents[1]] = [0 blockEvents[3:end]...]
            end

            if version_combined >= 1004000
                blockDurations[blockEvents[1]] = blockEvents[2]*blockDurationRaster
            else
                delayIDs_tmp[blockEvents[1]] = blockEvents[2]
            end
        end
    end
    sort(eventTable), sort(blockDurations), sort(delayIDs_tmp)
end

"""
read_blocks Read the [BLOCKS] section of a sequence file.
   library=read_blocks(fid) Read blocks from file identifier of an
   open MR sequence file and returns blockIDs, durations and delayIDs_tmp
   as arrays or vectors, rather than dicts

"""
function read_blocks_and_durs_as_arrays(io, blockDurationRaster, version_combined)
    if version_combined <= 1002001
        NumberBlockEvents = 7
    else
        NumberBlockEvents = 8
    end
    # we'll collect everything into a vector and then reshape later
    # we assume that we have at least 1000 blocks and allocate for that
    blocks = empty!(Vector{Int}(undef, NumberBlockEvents * 1000))  # go through empty! to preallocate
    blockDurations = empty!(Vector{Float64}(undef, 1000))  # go through empty! to preallocate
    delayIDs_tmp = empty!(Vector{Float64}(undef, 1000))  # go through empty! to preallocate
    num_lines = 0
    while true
        blockEvents = [parse(Int, d) for d in eachsplit(readline(io))]
        # todo: maybe better read chunks!? but then when do we stop? Maybe we should read entire file to ram before anyway...
        if length(blockEvents) == 0
            # empty line. break (fixme: pulseq doesn't forbid that, I believe, but other funcs here rely on that too for detecting end of blocks block...)
            # or maybe it forbids. p9 says 'subsequent lines...'
            break
        end
        if length(blockEvents) != NumberBlockEvents
            error("Expected $NumberBlockEvents events but got $(length(blockEvents)).")
        end
        if blockEvents[1] != 0  # id is not 0. Only then we care!
            if version_combined <= 1002001
                # add an extra 0 at the end to make equal to other version
                append!(blocks, blockEvents[3:end], 0) # only from 3 to end (discard id and duration, duration we treat below)
            else
                append!(blocks, blockEvents[3:end]) # only from 3 to end (discard id and duration, duration we treat below)
            end

            if version_combined >= 1004000
                push!(blockDurations, blockEvents[2] * blockDurationRaster)
            else
                push!(delayIDs_tmp, blockEvents[2])
            end
        end
        num_lines += 1
    end
    reshaped = reshape(blocks, :, num_lines)
    # we need 6 vals because we drop IDs and durations. earlier versions we added a 0 -> for all it's 8 - 2 = 6
    @assert size(reshaped, 1) == 6 "unexpected number of fields per block"  # for all versions!
    reshaped, blockDurations, delayIDs_tmp
end

"""
read_events Read an event section of a sequence file.
   library=read_events(fid) Read event data from file identifier of
   an open MR sequence file and return a library of events.

   library=read_events(fid,scale) Read event data and scale
   elements according to column vector scale.

   library=read_events(fid,scale,type) Attach the type string to
   elements of the library.

   library=read_events(...,library) Append new events to the given
   library.
"""
function read_events(io, scale; type=-1, eventLibrary=Dict())
    while true
        EventLength = length(scale) + 1
        fmt = Scanf.Format("%i"*"%f "^(EventLength-1))
        r, data... = scanf(readline(io), fmt, Int, zeros(Float64,EventLength-1)...)
        id = floor(Int, data[1])
        data = scale .* [data[2:end]...]'
        if type != -1
            eventLibrary[id] = Dict("data"=>data, "type"=>type)
        else
            eventLibrary[id] = Dict("data"=>data)
        end
        r == EventLength || break #Break on white space
    end
    sort(eventLibrary)
end

"""
read_shapes Read the [SHAPES] section of a sequence file.
   library=read_shapes(fid) Read shapes from file identifier of an
   open MR sequence file and return a library of shapes.
"""
function read_shapes(io, forceConvertUncompressed)
    shapeLibrary = Dict()
    readline(io)
    while true #Reading shapes
        r, id = @scanf(readline(io), "shape_id %i", Int)
        r == 1 || break #Break if no "shape_id id" is identified
        _, num_samples = @scanf(readline(io), "num_samples %i", Int)
        shape = Float64[]
        while true #Reading shape data
            r, data_point = @scanf(readline(io), "%f", Float64)
            r == 1 || break #Break if no sample
            append!(shape, data_point)
        end
        # check if conversion is needed: in v1.4.x we use length(data)==num_samples
        # as a marker for the uncompressed (stored) data. In older versions this condition could occur by chance
        if forceConvertUncompressed && length(shape)==num_samples
            shape = compress_shape(decompress_shape(shape, num_samples; forceDecompression=true))
            data = (num_samples, shape)
        else
            data = (num_samples, shape)
        end
        shapeLibrary[id] = data
    end

    sort(shapeLibrary)
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
        if forceCompression || num_samples > length(v)
            data = v
        else
            data = w
        end
    end

    num_samples, data
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
        for i = 1:length(dataPackMarkers)
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

    w
end

#"""
#READ Load sequence from file.
#   READ(seqObj, filename, ...) Read the given filename and load sequence
#   data into sequence object.
#
#   optional parwameter 'detectRFuse' can be given to let the function
#   infer the currently missing flags concerning the intended use of the RF
#   pulses (excitation, refocusing, etc). These are important for the
#   k-space trajectory calculation
#
#   Examples:
#   Load the sequence defined in gre.seq in my_sequences directory
#
#       read(seqObj,'my_sequences/gre.seq')
#
# See also  write
#"""
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
    version_combined = 0
    version_major = 0
    version_minor = 0
    gradLibrary = Dict()
    def = Dict()
    blockEvents = Dict()
    blockDurations = Dict()
    delayInd_tmp = Dict()
    rfLibrary = Dict()
    adcLibrary = Dict()
    tmp_delayLibrary = Dict()
    shapeLibrary = Dict()
    extensionLibrary = Dict()
    triggerLibrary = Dict()
    #Reading file and storing data
    open(filename) do io
        while !eof(io)
            section = readline(io)
            if      section == "[DEFINITIONS]"
                def = read_definitions(io)
            elseif  section == "[VERSION]"
                version_major, version_minor, _, version_combined = read_version(io)
            elseif  section == "[BLOCKS]"
                if version_combined == 0
                    @error "Pulseq file MUST include [VERSION] section prior to [BLOCKS] section"
                end
                blockEvents, blockDurations, delayInd_tmp = read_blocks(io, def["BlockDurationRaster"], version_combined)
            elseif  section == "[RF]"
                if version_combined >= 1004000
                    rfLibrary = read_events(io, [1/γ 1 1 1 1e-6 1 1]) # this is 1.4.x format
                else
                    rfLibrary = read_events(io, [1/γ 1 1 1e-6 1 1]) # this is 1.3.x and below
                    # we will have to scan through the library later after all the shapes have been loaded
                end
            elseif  section == "[GRADIENTS]"
                if version_combined >= 1004000
                    gradLibrary = read_events(io, [1/γ 1 1 1e-6]; type='g', eventLibrary=gradLibrary) # this is 1.4.x format
                else
                    gradLibrary = read_events(io, [1/γ 1 1e-6];   type='g', eventLibrary=gradLibrary) # this is 1.3.x and below
                end
            elseif  section == "[TRAP]"
                gradLibrary = read_events(io, [1/γ 1e-6 1e-6 1e-6 1e-6]; type='t', eventLibrary=gradLibrary);
            elseif  section == "[ADC]"
                adcLibrary = read_events(io, [1 1e-9 1e-6 1 1])
            elseif  section == "[DELAYS]"
                if version_combined >= 1004000
                    @error "Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section"
                end
                tmp_delayLibrary = read_events(io, 1e-6);
            elseif  section == "[SHAPES]"
                shapeLibrary = read_shapes(io, (version_major==1 && version_minor<4))
            elseif  section == "[EXTENSIONS]"
                extensionLibrary = read_events(io,[1 1 1]) #For now, it reads the extensions but it does not take it them into account
            elseif  section == "extension TRIGGERS 1"
                triggerLibrary = read_events(io,[1 1 1e-6 1e-6])
            elseif  section == "[SIGNATURE]"
                #Not implemented yet
            end

        end
    end
    # fix blocks, gradients and RF objects imported from older versions
    if version_combined < 1004000
        # scan through the RF objects
        for i = 0:length(rfLibrary)-1
            rfLibrary[i]["data"] = [rfLibrary[i]["data"][1:3]' 0 rfLibrary[i]["data"][4:end]']
        end
        # scan through the gradient objects and update 't'-s (trapezoids) und 'g'-s (free-shape gradients)
        for i = 0:length(gradLibrary)-1
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
            if gradLibrary[i]["type"] == 'g'
                #(1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
                gradLibrary[i]["data"] = [gradLibrary[i]["data"][1:2]; 0; gradLibrary[i]["data"][3:end]]
            end
        end
        # for versions prior to 1.4.0 blockDurations have not been initialized
        blockDurations = zeros(length(blockEvents))
        for i = 1:length(blockEvents)
        idelay = delayInd_tmp[i]
            if idelay > 0
                delay = tmp_delayLibrary[idelay]["data"][1]
                blockDurations[i] = delay
            end
        end
    end
    #Sequence
    obj = Dict(
        "blockEvents"=>blockEvents,
        "blockDurations"=>blockDurations,
        "delayInd_tmp"=>delayInd_tmp,
        "gradLibrary"=>gradLibrary,
        "rfLibrary"=>rfLibrary,
        "adcLibrary"=>adcLibrary,
        "tmp_delayLibrary"=>tmp_delayLibrary,
        "shapeLibrary"=>shapeLibrary,
        "extensionLibrary"=>extensionLibrary,
        "triggerLibrary"=>triggerLibrary,
        "definitions"=>def)
    #Transforming Dictionary to Sequence object
    #This should only work for Pulseq files >=1.4.0
    seq = Sequence()
    for i = 1:length(blockEvents)
        seq += get_block(obj,i)
    end
    # Final details
    # Remove dummy seq block at the start, Issue #203
    seq = seq[2:end]
    # Hack for including extension and triggers
    seq.DEF["additional_text"] = read_Extension(extensionLibrary, triggerLibrary) #Temporary hack
    seq.DEF = recursive_merge(obj["definitions"], seq.DEF)
    # Koma specific details for reconstrucion
    seq.DEF["FileName"] = basename(filename)
    seq.DEF["PulseqVersion"] = version_combined
    if !haskey(seq.DEF,"Nx")
        Nx = maximum(seq.ADC.N)
        RF_ex = (get_flip_angles(seq) .<= 90.01) .* is_RF_on.(seq)
        Nz = max(length(unique(seq.RF[RF_ex].Δf)), 1)
        Ny = sum(is_ADC_on.(seq)) / Nz |> x->floor(Int,x)

        seq.DEF["Nx"] = Nx  #Number of samples per ADC
        seq.DEF["Ny"] = Ny  #Number of ADC events
        seq.DEF["Nz"] = Nz  #Number of unique RF frequencies, in a 3D acquisition this should not work
    end
    #Koma sequence
    return seq
end

#To Sequence
"""
    grad = read_Grad(gradLibrary, shapeLibrary, Δt_gr, i)

Reads the gradient. It is used internally by [`get_block`](@ref).

# Arguments
- `gradLibrary`: (`::Dict{K, V}`) the "gradLibrary" dictionary
- `shapeLibrary`: (`::Dict{K, V}`) the "shapeLibrary" dictionary
- `Δt_gr`: (`::Float64`, `[s]`) gradient raster time
- `i`: (`::Int64`) index of the axis in the block event

# Returns
- `grad`: (::Grad) Gradient struct
"""
function read_Grad(gradLibrary, shapeLibrary, Δt_gr, i)
    G = Grad(0,0)
    if gradLibrary[i]["type"] == 't' #if trapezoidal gradient
        #(1)amplitude (2)rise (3)flat (4)fall (5)delay
        g_A, g_rise, g_T, g_fall, g_delay = gradLibrary[i]["data"]
        G = Grad(g_A,g_T,g_rise,g_fall,g_delay)
    elseif gradLibrary[i]["type"] == 'g' #Arbitrary gradient waveform
        #(1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
        g = gradLibrary[i]["data"]
        amplitude =     g[1]
        amp_shape_id =  g[2] |> x->floor(Int64,x)
        time_shape_id = g[3] |> x->floor(Int64,x)
        delay =         g[4]
        #Amplitude
        gA = amplitude * decompress_shape(shapeLibrary[amp_shape_id]...)
        Nrf = length(gA) - 1
        #Creating timings
        if time_shape_id == 0 #no time waveform
            gT = Nrf * Δt_gr
            G = Grad(gA, gT, Δt_gr/2, Δt_gr/2, delay)
        else
            gt = decompress_shape(shapeLibrary[time_shape_id]...)
            gT = (gt[2:end] .- gt[1:end-1]) * Δt_gr
            G = Grad(gA,gT,0,0,delay)
        end
    end
    G
end

"""
    rf = read_RF(rfLibrary, shapeLibrary, Δt_rf, i)

Reads the RF. It is used internally by [`get_block`](@ref).

# Arguments
- `rfLibrary`: (`::Dict{K, V}`) the "rfLibrary" dictionary
- `shapeLibrary`: (`::Dict{K, V}`) the "shapeLibrary" dictionary
- `Δt_rf`: (`::Float64`, `[s]`) RF raster time
- `i`: (`::Int64`) index of the RF in the block event

# Returns
- `rf`: (`1x1 ::Matrix{RF}`) RF struct
"""
function read_RF(rfLibrary, shapeLibrary, Δt_rf, i)
    #Unpacking
    #(1)amplitude (2)mag_id (3)phase_id (4)time_shape_id (5)delay (6)freq (7)phase
    r = rfLibrary[i]["data"]
    amplitude =     r[1]
    mag_id =        r[2] |> x->floor(Int64,x)
    phase_id =      r[3] |> x->floor(Int64,x)
    time_shape_id = r[4] |> x->floor(Int64,x)
    delay =         r[5] + (time_shape_id==0)*Δt_rf/2
    freq =          r[6]
    phase =         r[7]
    #Amplitude and phase waveforms
    if amplitude != 0 && mag_id != 0
        rfA = decompress_shape(shapeLibrary[mag_id]...)[1:end-1]
        rfϕ = decompress_shape(shapeLibrary[phase_id]...)[1:end-1]
        @assert all(rfϕ.>=0) "[RF id $i] Phase waveform rfϕ must have non-negative samples (1.>=rfϕ.>=0). "
        Nrf = shapeLibrary[mag_id][1] - 1
        rfAϕ = amplitude .* rfA .* exp.(-1im*(2π*rfϕ .+ phase))
    else
        rfA = 0
        rfϕ = 0
        Nrf = 0
        rfAϕ = 0
    end
    #Creating timings
    if time_shape_id == 0 #no time waveform
        rfT = Nrf * Δt_rf
    else
        rft = decompress_shape(shapeLibrary[time_shape_id]...)
        rfT = (rft[2:end] .- rft[1:end-1]) * Δt_rf
    end
    R = reshape([RF(rfAϕ,rfT,freq,delay)],1,1)#[RF(rfAϕ,rfT,freq,delay);;]
    R
end

"""
    adc = read_ADC(adcLibrary, i)

Reads the ADC. It is used internally by [`get_block`](@ref).

# Arguments
- `adcLibrary`: (`::Dict{String, Any}`) the "adcLibrary" dictionary
- `i`: (`::Int64`) index of the adc in the block event

# Returns
- `adc`: (`1x1 ::Vector{ADC}`) ADC struct
"""
function read_ADC(adcLibrary, i)
    #Unpacking
    #(1)num (2)dwell (3)delay (4)freq (5)phase
    if !isempty(adcLibrary) # Is this the best? maybe defining i=0 is better, it works with RFs(?)
        a = adcLibrary[i]["data"]
    else
        a = [0,0,0,0,0]
    end
    num =   a[1] |> x->floor(Int64,x)
    dwell = a[2]
    delay = a[3] + dwell/2
    freq =  a[4]
    phase = a[5]
    #Definition
    T = (num-1) * dwell
    A = [ADC(num,T,delay,freq,phase)]
    A
end

"""
    ext = read_Extension(extensionLibrary, triggerLibrary, i)

Reads the Extension. It is used internally by [`get_block`](@ref).

# Arguments
- `extensionLibrary`: (`::Dict{String, Any}`) the "extensionLibrary" dictionary
- `i`: (`::Int64`) index of the ext in the block event

# Returns
- `ext`: (`1x1 ::Vector{ADC}`) Extension struct
"""
function read_Extension(extensionLibrary, triggerLibrary)
    # Only uses triggers and one extension per block
    # Unpacking
    # Extensions
    # (1)id (2)type (3)ref (4)next_id
    # if !isempty(extensionLibrary)
    #     e = extensionLibrary[i]["data"]
    # else
    #     e = [0,0,0,0]
    # end
    # type    = e[1] |> x->floor(Int64,x) #1=Trigger
    # ref     = e[2] |> x->floor(Int64,x)
    # next_id = e[3] |> x->floor(Int64,x)
    # Trigger
    # (1)id (2)type (3)channel (4)delay (5)duration
    # if !isempty(triggerLibrary)
    #     t = triggerLibrary[ref]["data"]
    # else
    #     t = [0,0,0,0]
    # end
    # type     = t[1] |> x->floor(Int64,x)
    # channel  = t[2] |> x->floor(Int64,x)
    # delay    = t[3]
    # duration = t[4]
    #Definition
    # trig = Trigger(type, channel, delay, duration)
    # E = [Extension([trig])]
    # E = Dict("extension"=>[ref]) 
    additional_text =  
    """# Format of extension lists:
    # id type ref next_id
    # next_id of 0 terminates the list
    # Extension list is followed by extension specifications
    [EXTENSIONS]
    """
    for id = eachindex(extensionLibrary)
        (id == 0) && continue 
        type, ref, next_id = floor.(Int64, extensionLibrary[id]["data"])
        additional_text *= "$id $type $ref $next_id\n" 
    end
    additional_text *=
    """

    # Extension specification for digital output and input triggers:
    # id type channel delay (us) duration (us)
    extension TRIGGERS 1
    """
    for id = eachindex(triggerLibrary)
        (id == 0) && continue 
        type, channel, delay, duration = floor.(Int64, round.([1 1 1e6 1e6] .* triggerLibrary[id]["data"]))
        additional_text *= "$id $type $channel $delay $duration\n" 
    end
    return additional_text
end

"""
    seq = get_block(obj, i)

Block sequence definition. Used internally by [`read_seq`](@ref).

# Arguments
- `obj`: (`::Dict{String, Any}`) main dictionary
- `i`: (`::Int64`) index of a block event

# Returns
- `s`: (`::Sequence`) block Sequence struct
"""
function get_block(obj, i)
    #Unpacking
    idelay, irf, ix, iy, iz, iadc, iext = obj["blockEvents"][i]

    #Gradient definition
    Δt_gr = obj["definitions"]["GradientRasterTime"]
    Gx = read_Grad(obj["gradLibrary"], obj["shapeLibrary"], Δt_gr, ix)
    Gy = read_Grad(obj["gradLibrary"], obj["shapeLibrary"], Δt_gr, iy)
    Gz = read_Grad(obj["gradLibrary"], obj["shapeLibrary"], Δt_gr, iz)
    G = reshape([Gx;Gy;Gz],3,1) #[Gx;Gy;Gz;;]

    #RF definition
    Δt_rf = obj["definitions"]["RadiofrequencyRasterTime"]
    R = read_RF(obj["rfLibrary"], obj["shapeLibrary"], Δt_rf, irf)

    #ADC definition
    A = read_ADC(obj["adcLibrary"], iadc)

    #DUR
    D = [obj["blockDurations"][i]]

    #Extensions
    E = Dict("extension"=>[iext]) #read_Extension(obj["extensionLibrary"], iext, i)

    #Sequence block definition
    s = Sequence(G,R,A,D,E)
    s
end

"""
does what get_blocks does, but for many blocks at once
and works on array that contains block ids directly
it will use all blocks (2nd dim in blockEvents) to create the resulting sequence
"""
function get_seq_from_blocks(libraries, blockEvents::Array{Int, 2}, blockDurations::Vector{Float64})
    # get some views for blocks things
    ids_rf = blockEvents[1, :]
    ids_gradx = blockEvents[2, :]
    ids_grady = blockEvents[3, :]
    ids_gradz = blockEvents[4, :]
    ids_adc = blockEvents[5, :]
    ids_ext = blockEvents[6, :]

    # grads first
    Δt_gr = libraries["definitions"]["GradientRasterTime"]
    # allocate grads array space
    num_blocks = size(blockEvents, 2)
    GR = Array{Grad}(undef, 3, num_blocks)
    GR[1, :] = collect(read_Grad(libraries["gradLibrary"], libraries["shapeLibrary"], Δt_gr, id) for id in ids_gradx)
    GR[2, :] = collect(read_Grad(libraries["gradLibrary"], libraries["shapeLibrary"], Δt_gr, id) for id in ids_grady)
    GR[3, :] = collect(read_Grad(libraries["gradLibrary"], libraries["shapeLibrary"], Δt_gr, id) for id in ids_gradz)

    #RFs
    Δt_rf = libraries["definitions"]["RadiofrequencyRasterTime"]
    # result of read_RF is a 1-element matrix...
    RFs = empty!(Vector{RF}(undef, num_blocks))
    append!(RFs, read_RF(libraries["rfLibrary"], libraries["shapeLibrary"], Δt_rf, id)[1, 1] for id in ids_rf)
    RFs = reshape(RFs, 1, num_blocks)

    # ADC
    ADCs = empty!(Vector{ADC}(undef, num_blocks))
    # result of adc is a 1-element vector...
    append!(ADCs, (read_ADC(libraries["adcLibrary"], id)[1] for id in ids_adc))

    #DEFs which here are only Extensions. we can collect them in a vector
    # and then set the dict key
    DEFs = Dict("extension"=>ids_ext)
    Sequence(GR, RFs, ADCs, blockDurations, DEFs)
end

"""

    seq = read_seq_via_blocks_as_int_array(filename)

Returns the Sequence struct from a Pulseq file with `.seq` extension.

ALTERNATIVE IMPLEMENTATION OF read_seq.
This reads blocks as a big int array. Does not sort blocks by ID!
(pulseq specification says not to order IDs, although it's done in matlab!)
Calls costly Sequence constructor only once.

# Arguments
- `filename`: (`::String`) absolute or relative path of the sequence file `.seq`

# Returns
- `seq`: (`::Sequence`) Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq_via_blocks_as_int_array(seq_file)

julia> plot_seq(seq)
```
"""
function read_seq_via_blocks_as_int_array(filename)
    @info "Loading sequence $(basename(filename)) ..."
    version_combined = 0
    version_major = 0
    version_minor = 0
    gradLibrary = Dict()
    def = Dict()
    rfLibrary = Dict()
    adcLibrary = Dict()
    tmp_delayLibrary = Dict()
    shapeLibrary = Dict()
    extensionLibrary = Dict()
    triggerLibrary = Dict()
    blockEvents = Array{Int, 2}(undef, 0, 0)
    blockDurations = Vector{Float64}(undef, 0)
    delayInd_tmpVec = Vector{Int}(undef, 0)
    #Reading file and storing data
    open(filename) do io
        while !eof(io)
            section = readline(io)
            if      section == "[DEFINITIONS]"
                def = read_definitions(io)
            elseif  section == "[VERSION]"
                version_major, version_minor, _, version_combined = read_version(io)
            elseif  section == "[BLOCKS]"
                if version_combined == 0
                    @error "Pulseq file MUST include [VERSION] section prior to [BLOCKS] section"
                end
                blockEvents, blockDurations, delayInd_tmpVec = read_blocks_and_durs_as_arrays(io, def["BlockDurationRaster"], version_combined)
            elseif  section == "[RF]"
                if version_combined >= 1004000
                    rfLibrary = read_events(io, [1/γ 1 1 1 1e-6 1 1]) # this is 1.4.x format
                else
                    rfLibrary = read_events(io, [1/γ 1 1 1e-6 1 1]) # this is 1.3.x and below
                    # we will have to scan through the library later after all the shapes have been loaded
                end
            elseif  section == "[GRADIENTS]"
                if version_combined >= 1004000
                    gradLibrary = read_events(io, [1/γ 1 1 1e-6]; type='g', eventLibrary=gradLibrary) # this is 1.4.x format
                else
                    gradLibrary = read_events(io, [1/γ 1 1e-6];   type='g', eventLibrary=gradLibrary) # this is 1.3.x and below
                end
            elseif  section == "[TRAP]"
                gradLibrary = read_events(io, [1/γ 1e-6 1e-6 1e-6 1e-6]; type='t', eventLibrary=gradLibrary);
            elseif  section == "[ADC]"
                adcLibrary = read_events(io, [1 1e-9 1e-6 1 1])
            elseif  section == "[DELAYS]"
                if version_combined >= 1004000
                    @error "Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section"
                end
                tmp_delayLibrary = read_events(io, 1e-6);
            elseif  section == "[SHAPES]"
                shapeLibrary = read_shapes(io, (version_major==1 && version_minor<4))
            elseif  section == "[EXTENSIONS]"
                extensionLibrary = read_events(io,[1 1 1]) #For now, it reads the extensions but it does not take it them into account
            elseif  section == "extension TRIGGERS 1"
                triggerLibrary = read_events(io,[1 1 1e-6 1e-6])
            elseif  section == "[SIGNATURE]"
                #Not implemented yet
            end

        end
    end
    # fix blocks, gradients and RF objects imported from older versions
    if version_combined < 1004000
        # scan through the RF objects
        for i = 0:length(rfLibrary)-1
            rfLibrary[i]["data"] = [rfLibrary[i]["data"][1:3]' 0 rfLibrary[i]["data"][4:end]']
        end
        # scan through the gradient objects and update 't'-s (trapezoids) und 'g'-s (free-shape gradients)
        for i = 0:length(gradLibrary)-1
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
            if gradLibrary[i]["type"] == 'g'
                #(1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
                gradLibrary[i]["data"] = [gradLibrary[i]["data"][1:2]; 0; gradLibrary[i]["data"][3:end]]
            end
        end
        # for versions prior to 1.4.0 blockDurations have not been initialized
        resize!(blockDurations, size(blockEvents, 2))
        for i = 1:size(blockEvents, 1)
        idelay = delayInd_tmpVec[i]
            if idelay > 0
                delay = tmp_delayLibrary[idelay]["data"][1]
                blockDurations[i] = delay
            end
        end
    end
    #Sequence
    libraries = Dict(
        "gradLibrary"=>gradLibrary,
        "rfLibrary"=>rfLibrary,
        "adcLibrary"=>adcLibrary,
        "tmp_delayLibrary"=>tmp_delayLibrary,
        "shapeLibrary"=>shapeLibrary,
        "extensionLibrary"=>extensionLibrary,
        "triggerLibrary"=>triggerLibrary,
        "definitions"=>def)
    #Transforming Dictionary to Sequence object
    #This should only work for Pulseq files >=1.4.0
    # in this most optimized implementation, the following lines still take
    # 45% of all the time. These have to be optimized, as we still need
    # 2s for a 2d 18point mrf seq with 64 PE lines.
    # a 3d seq with 300 points and say 30x30 = 900 lines 
    # would take 234x longer -> 8min!

    seq = get_seq_from_blocks(libraries, blockEvents, blockDurations)
    # old code
    # seq = Sequence()
    # for i = 1:length(blockEvents)
    #     seq += get_block(obj,i)
    # end
    # Final details
    # Remove dummy seq block at the start, Issue #203
    # not needed anymore
    # seq = seq[2:end]
    # Hack for including extension and triggers
    seq.DEF["additional_text"] = read_Extension(extensionLibrary, triggerLibrary) #Temporary hack
    seq.DEF = recursive_merge(libraries["definitions"], seq.DEF)
    # Koma specific details for reconstrucion
    seq.DEF["FileName"] = basename(filename)
    seq.DEF["PulseqVersion"] = version_combined
    if !haskey(seq.DEF,"Nx")
        Nx = maximum(seq.ADC.N)
        RF_ex = (get_flip_angles(seq) .<= 90.01) .* is_RF_on.(seq)
        Nz = max(length(unique(seq.RF[RF_ex].Δf)), 1)
        Ny = sum(is_ADC_on.(seq)) / Nz |> x->floor(Int,x)

        seq.DEF["Nx"] = Nx  #Number of samples per ADC
        seq.DEF["Ny"] = Ny  #Number of ADC events
        seq.DEF["Nz"] = Nz  #Number of unique RF frequencies, in a 3D acquisition this should not work
    end
    #Koma sequence
    return seq
end
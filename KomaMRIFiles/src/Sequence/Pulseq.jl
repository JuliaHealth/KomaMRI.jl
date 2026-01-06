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
        @warn "Pulseq $(pulseq_version) not yet supported by this KomaMRIFiles release. Track progress at https://github.com/JuliaHealth/KomaMRI.jl/pull/614"
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
    if !haskey(def,"BlockDurationRaster")       def["BlockDurationRaster"] = 1e-5       end
    if !haskey(def,"GradientRasterTime")        def["GradientRasterTime"] = 1e-5        end
    if !haskey(def,"RadiofrequencyRasterTime")  def["RadiofrequencyRasterTime"] = 1e-6  end
    if !haskey(def,"AdcRasterTime")             def["AdcRasterTime"] = 1e-7             end
    def
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
    return isnothing(signature_type) && isnothing(signature_hash) ? nothing : (type = signature_type, hash = signature_hash)
end

"""
read_blocks Read the [BLOCKS] section of a sequence file.
   library=read_blocks(fid) Read blocks from file identifier of an
   open MR sequence file and return the event table.
"""
function read_blocks(io, blockDurationRaster, pulseq_version)
    eventTable = Dict{Int64, Vector{Int64}}()
    blockDurations = Dict{Int64, Float64}()
    delayIDs_tmp = Dict{Int64, Float64}()
    while true
        if pulseq_version <= v"1.2.1"
            NumberBlockEvents = 7
        else
            NumberBlockEvents = 8
        end

        read_event = readline(io)
        !isempty(read_event) || break
        blockEvents = parse.(Int64, split(read_event))

        if blockEvents[1] != 0
            if pulseq_version <= v"1.2.1"
                eventTable[blockEvents[1]] = Int64[0; blockEvents[3:end]...; 0]
            else
                eventTable[blockEvents[1]] = Int64[0; blockEvents[3:end]...]
            end

            if pulseq_version >= v"1.4.0"
                blockDurations[blockEvents[1]] = blockEvents[2]*blockDurationRaster
            else
                delayIDs_tmp[blockEvents[1]] = blockEvents[2]
            end
        end

        length(blockEvents) == NumberBlockEvents || break #Break on white space
    end
    eventTable, blockDurations, delayIDs_tmp
end

"""
read_events Read an event section of a sequence file.
   library=read_events(fid,scale,format) Read event data and scale
   elements according to column vector scale and format string.

   library=read_events(...; eventLibrary=library) Append new events to the given
   library.
"""
function read_events(io, scale; format="%i "*"%f "^(length(scale)), eventLibrary=Dict())
    eventLength = length(scale) + 1
    fmt = Scanf.Format(format)
    args = Tuple([f == "%i" ? Int : (f == "%f" ? Float64 : (f == "%c" ? Char : String)) for f in split(format)])
    while true
        line = readline(io)
        isempty(line) && break
        r, data... = scanf(line, fmt, args...)
        r == eventLength || break #Break if not all values read
        id = floor(Int, data[1])
        data = [d isa Char ? d : s*d for (s, d) in zip(scale, data[2:end])]   
        eventLibrary[id] = Dict("data"=>data)
    end
    eventLibrary
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
            data_point = tryparse(Float64, readline(io))
            !isnothing(data_point) || break #Break if no sample
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
    shapeLibrary
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

    w
end

"""
    fix_first_last_grads!(seq::Sequence)

Updates the Sequence `seq` with new first and last points for gradients.
"""
function fix_first_last_grads!(seq::Sequence)
    # Add first and last Pulseq points
    grad_prev_last = [0.0; 0.0; 0.0]
    for bi in 1:length(seq)
        for gi in 1:3
            gr = seq.GR[gi, bi]
            if seq.DUR[bi] > 0
                if sum(abs.(gr.A)) == 0   # this is for no-gradient case
                    grad_prev_last[gi] = 0.0
                    continue
                else
                    if gr.A isa Vector # this is for the uniformly-shaped or time-shaped case
                        if gr.delay > 0
                            grad_prev_last[gi] = 0.0
                        end
                        seq.GR[gi, bi].first = grad_prev_last[gi]
                        if gr.T isa Array # this is for time-shaped case (I assume it is the extended trapezoid case)
                            seq.GR[gi, bi].last = gr.A[end]
                        else
                            odd_step1 = [seq.GR[gi, bi].first; 2 * gr.A]
                            odd_step2 = odd_step1 .* (mod.(1:length(odd_step1), 2) * 2 .- 1)
                            waveform_odd_rest = cumsum(odd_step2) .* (mod.(1:length(odd_step2), 2) * 2 .- 1)
                            seq.GR[gi, bi].last = waveform_odd_rest[end]
                        end
                        grad_prev_last[gi] = seq.GR[gi, bi].last
                    else    # this is for the trapedoid case
                        grad_prev_last[gi] = 0.0
                    end
                end
            end
        end
    end
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
    gradLibrary = Dict()
    trapLibrary = Dict()
    def = Dict()
    signature = nothing
    blockEvents = Dict()
    blockDurations = Dict()
    delayInd_tmp = Dict()
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
                blockEvents, blockDurations, delayInd_tmp = read_blocks(io, def["BlockDurationRaster"], pulseq_version)
            elseif  section == "[RF]"
                if pulseq_version >= v"1.5.0"
                    rfLibrary = read_events(io, [1/γ 1 1 1 1 1e-6 1 1 1 1 1]; format="%i "*"%f "^(10)*"%c ") # this is 1.5.x format
                elseif pulseq_version >= v"1.4.0"
                    rfLibrary = read_events(io, [1/γ 1 1 1 1e-6 1 1]) # this is 1.4.x format
                else
                    rfLibrary = read_events(io, [1/γ 1 1 1e-6 1 1]) # this is 1.3.x and below
                    # we will have to scan through the library later after all the shapes have been loaded
                end
            elseif  section == "[GRADIENTS]"
                if pulseq_version >= v"1.5.0"
                    gradLibrary = read_events(io, [1/γ 1/γ 1/γ 1 1 1e-6]; eventLibrary=gradLibrary) # this is 1.5.x format
                elseif pulseq_version >= v"1.4.0"
                    gradLibrary = read_events(io, [1/γ 1 1 1e-6]; eventLibrary=gradLibrary) # this is 1.4.x format
                else
                    gradLibrary = read_events(io, [1/γ 1 1e-6]; eventLibrary=gradLibrary) # this is 1.3.x and below
                end
            elseif  section == "[TRAP]"
                trapLibrary = read_events(io, [1/γ 1e-6 1e-6 1e-6 1e-6]; eventLibrary=trapLibrary);
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
            else
                if startswith(section, "extension")
                    ext = section[11:end]
                    ext_string = split(ext, " ")[1]
                    ext_type   = KomaMRIBase.get_EXT_type_from_symbol(Val(Symbol(ext_string)))
                    ext_id     = parse(Int, split(ext, " ")[2])
                    if !haskey(def,"RequiredExtensions") def["RequiredExtensions"] = []  end

                    if ext_type <: Extension
                        extensionTypeLibrary[ext_id] = ext_type
                        extensionSpecLibrary[ext_id] = read_events(io, KomaMRIBase.get_scale(ext_type); format="%i "*KomaMRIBase.get_scanf_format(ext_type))
                    else
                        if ext_string in def["RequiredExtensions"]
                            @error "Extension $ext_string is required by the sequence but not supported by KomaMRI reader"
                        else
                            @warn "Ignoring unsupported extension: $ext_string"
                        end
                    end
                else
                    @error "Unknown section code: $section"
                end
            end
        end
    end
    verify_signature!(filename, signature; pulseq_version=pulseq_version)
    # fix blocks, gradients and RF objects imported from older versions
    if pulseq_version < v"1.4.0"
        # scan through the RF objects
        for i in eachindex(rfLibrary)
            rfLibrary[i]["data"] = [rfLibrary[i]["data"][1:3]; 0.0; rfLibrary[i]["data"][4:end]]
        end
        # scan through the trapezoidal gradient objects
        for i in eachindex(trapLibrary)
            #(1)amplitude (2)rise (2)flat (3)fall (4)delay
            if trapLibrary[i]["data"][2] == 0 #rise
                if abs(trapLibrary[i]["data"][1]) == 0 && trapLibrary[i]["data"][3] > 0
                    trapLibrary[i]["data"][3] -= def["gradRasterTime"]
                    trapLibrary[i]["data"][2]  = def["gradRasterTime"]
                end
            end
            if trapLibrary[i]["data"][4] == 0 #delay
                if abs(trapLibrary[i]["data"][1]) == 0 && trapLibrary[i]["data"][3] > 0
                    trapLibrary[i]["data"][3] -= def["gradRasterTime"]
                    trapLibrary[i]["data"][4]  = def["gradRasterTime"]
                end
            end
        end
        # scan through the time-shaped gradient objects
        for i in eachindex(gradLibrary)
            #(1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
            gradLibrary[i]["data"] = [gradLibrary[i]["data"][1:2]; 0.0; gradLibrary[i]["data"][3:end]]
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
        "trapLibrary"=>trapLibrary,
        "rfLibrary"=>rfLibrary,
        "adcLibrary"=>adcLibrary,
        "tmp_delayLibrary"=>tmp_delayLibrary,
        "shapeLibrary"=>shapeLibrary,
        "extensionInstanceLibrary"=>extensionInstanceLibrary,
        "extensionTypeLibrary"=>extensionTypeLibrary,
        "extensionSpecLibrary"=>extensionSpecLibrary,
        "definitions"=>def)
    #Transforming Dictionary to Sequence object
    #This should only work for Pulseq files >=1.4.0
    seq = Sequence()
    for i = 1:length(blockEvents)
        seq += get_block(obj, i)
    end
    # Add first and last points for gradients #320 for version <= 1.4.2
    if pulseq_version < v"1.5.0"
        fix_first_last_grads!(seq)
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
    grad = read_Grad(gradLibrary, trapLibrary, shapeLibrary, Δt_gr, i)

Reads the gradient. It is used internally by [`get_block`](@ref).

# Arguments
- `gradLibrary`: (`::Dict{K, V}`) the "gradLibrary" dictionary
- `trapLibrary`: (`::Dict{K, V}`) the "trapLibrary" dictionary
- `shapeLibrary`: (`::Dict{K, V}`) the "shapeLibrary" dictionary
- `Δt_gr`: (`::Float64`, `[s]`) gradient raster time
- `i`: (`::Int64`) index of the axis in the block event

# Returns
- `grad`: (::Grad) Gradient struct
"""
function read_Grad(gradLibrary, trapLibrary, shapeLibrary, Δt_gr, i)
    G = Grad(0.0,0.0)
    if haskey(trapLibrary, i) # Trapezoidal gradient waveform
        #(1)amplitude (2)rise (3)flat (4)fall (5)delay
        g_A, g_rise, g_T, g_fall, g_delay = trapLibrary[i]["data"]
        G = Grad(g_A,g_T,g_rise,g_fall,g_delay,0.0,0.0)
    elseif haskey(gradLibrary, i) # Arbitrary gradient waveform
        g = gradLibrary[i]["data"]
        if length(g) == 6 # for version 1.5.x
            #(1)amplitude (2)first_grads (3)last_grads (4)amp_shape_id (5)time_shape_id (6)delay
            amplitude     = g[1]
            first_grads   = g[2]
            last_grads    = g[3]
            amp_shape_id  = g[4] |> x->floor(Int64,x)
            time_shape_id = g[5] |> x->floor(Int64,x)
            delay         = g[6]
        else # for version 1.4.x and below
            #(1)amplitude (2)amp_shape_id (3)time_shape_id (4)delay
            amplitude     = g[1]
            amp_shape_id  = g[2] |> x->floor(Int64,x)
            time_shape_id = g[3] |> x->floor(Int64,x)
            delay         = g[4]
            first_grads   = 0.0
            last_grads    = 0.0
        end
        #Amplitude
        gA = amplitude * decompress_shape(shapeLibrary[amp_shape_id]...)
        Ngr = length(gA) - 1
        #Creating timings
        if time_shape_id == 0 #no time waveform. Default time raster
            gT = Ngr * Δt_gr
            G = Grad(gA, gT, Δt_gr/2, Δt_gr/2, delay, first_grads, last_grads)
        elseif time_shape_id == -1 #New in pulseq 1.5.x: no time waveform. 1/2 of the default time raster
            gT = Ngr * Δt_gr / 2
            G = Grad(gA, gT, Δt_gr/2, Δt_gr/2, delay, first_grads, last_grads)
        else
            gt = decompress_shape(shapeLibrary[time_shape_id]...)
            delay += gt[1] * Δt_gr # offset due to the shape starting at a non-zero value
            gT = diff(gt) * Δt_gr
            G = Grad(gA,gT,0.0,0.0,delay,first_grads,last_grads)
        end
    end
    return G
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
    if isempty(rfLibrary) || i==0
        return reshape([RF(0.0,0.0)], 1, 1)
    end
    #Unpacking
    r = rfLibrary[i]["data"]
    if length(r) == 11 # for version 1.5.x
        #(1)amplitude (2)mag_id (3)phase_id (4)time_shape_id (5)center (6)delay (7)freq_ppm (8)phase_ppm (9)freq (10)phase (11)use
        amplitude     = r[1]
        mag_id        = r[2] |> x->floor(Int64,x)
        phase_id      = r[3] |> x->floor(Int64,x)
        time_shape_id = r[4] |> x->floor(Int64,x)
        center        = r[5]
        delay         = r[6] + (time_shape_id==0)*Δt_rf/2
        freq_ppm      = r[7]
        phase_ppm     = r[8]
        freq          = r[9]
        phase         = r[10]
        use           = r[11]
    else # for version 1.4.x and below
        #(1)amplitude (2)mag_id (3)phase_id (4)time_shape_id (5)delay (6)freq (7)phase
        amplitude     = r[1]
        mag_id        = r[2] |> x->floor(Int64,x)
        phase_id      = r[3] |> x->floor(Int64,x)
        time_shape_id = r[4] |> x->floor(Int64,x)
        delay         = r[5] + (time_shape_id==0)*Δt_rf/2
        freq          = r[6]
        phase         = r[7]
        center        = 0.0
        use           = 'u'
    end
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
    use = KomaMRIBase.get_RF_use_from_char(Val(use))
    if length(r) == 11 # for version 1.5.x
        return [RF(rfAϕ,rfT,freq,delay,center,use);;]
    else # for version 1.4.x and below
        return [RF(rfAϕ,rfT,freq,delay);;]
    end
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

    if isempty(adcLibrary) || i==0
        return [ADC(0, 0)]
    end

    #Unpacking
    a = adcLibrary[i]["data"]
    if length(a) == 8 # for version 1.5.x
        #(1)num (2)dwell (3)delay (4)freq_ppm (5)phase_ppm (6)freq (7)phase (8)phase_shape_id
        num       = a[1] |> x->floor(Int64,x)
        dwell     = a[2]
        delay     = a[3] + dwell/2
        freq_ppm  = a[4]
        phase_ppm = a[5]
        freq      = a[6]
        phase     = a[7]
        phase_id  = a[8] |> x->floor(Int64,x)
    else # for version 1.4.x and below
        #(1)num (2)dwell (3)delay (4)freq (5)phase
        num       = a[1] |> x->floor(Int64,x)
        dwell     = a[2]
        delay     = a[3] + dwell/2
        freq      = a[4]
        phase     = a[5]
        freq_ppm  = 0.0
        phase_ppm = 0.0
        phase_id  = 0
    end
    #Definition
    T = (num-1) * dwell
    return [ADC(num,T,delay,freq,phase)]
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
    idur, irf, igx, igy, igz, iadc, iext = obj["blockEvents"][i]

    #Gradient definition
    Δt_gr = obj["definitions"]["GradientRasterTime"]
    Gx = read_Grad(obj["gradLibrary"], obj["trapLibrary"], obj["shapeLibrary"], Δt_gr, igx)
    Gy = read_Grad(obj["gradLibrary"], obj["trapLibrary"], obj["shapeLibrary"], Δt_gr, igy)
    Gz = read_Grad(obj["gradLibrary"], obj["trapLibrary"], obj["shapeLibrary"], Δt_gr, igz)
    G = reshape([Gx;Gy;Gz],3,1) #[Gx;Gy;Gz;;]

    #RF definition
    Δt_rf = obj["definitions"]["RadiofrequencyRasterTime"]
    R = read_RF(obj["rfLibrary"], obj["shapeLibrary"], Δt_rf, irf)

    #ADC definition
    A = read_ADC(obj["adcLibrary"], iadc)

    #DUR
    D = Float64[max(obj["blockDurations"][i], dur(Gx), dur(Gy), dur(Gz), dur(R[1]), dur(A[1]))]

    #Extensions
    E = read_extension(obj["extensionInstanceLibrary"], obj["extensionTypeLibrary"], obj["extensionSpecLibrary"], iext)

    # Definitition
    DEF = Dict{String,Any}()

    #Sequence block definition
    s = Sequence(G,R,A,D,E,DEF)
    s
end

"""
    EXT = read_extension(extensionInstanceLibrary, extensionTypeLibrary, extensionSpecLibrary, i)

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
function read_extension(extensionInstanceLibrary, extensionTypeLibrary, extensionSpecLibrary, i)
    EXT = [Extension[]]
    if isempty(extensionInstanceLibrary) || i==0
        return EXT
    end
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
    lowercase(bytes2hex(digest))
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

function verify_signature!(filename::String, signature::Union{Nothing,NamedTuple}; pulseq_version::VersionNumber=v"1.4.0")
    isnothing(signature) && begin
        @warn "Pulseq [SIGNATURE] section is missing; skipping verification." filename
        return
    end
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
        # Don't error, just warn - the file can still be used
    end
end

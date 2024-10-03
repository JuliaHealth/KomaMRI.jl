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
read_definitions Read the [SIGNATURE] section of a sequence file.
   signature=read_signature(fid) Read user signature from file
   identifier of an open MR sequence file and return a map of
   key/value entries.
"""
function read_signature(io)
    signature = ""
    while true
        line = readline(io)
        line_split = String.(split(line))
        (length(line_split) > 0) || break #Break on white space
        key = line_split[1]
        if (key == "Hash")
            value_string_array = line_split[2:end]
            parsed_array = [tryparse(Float64, s) === nothing ? s : tryparse(Float64, s) for s = value_string_array]
            signature = length(parsed_array) == 1 ? parsed_array[1] : parsed_array
        end
    end
    return signature
end

"""
read_blocks Read the [BLOCKS] section of a sequence file.
   library=read_blocks(fid) Read blocks from file identifier of an
   open MR sequence file and return the event table.
"""
function read_blocks(io, blockDurationRaster, version_combined)
    eventTable = Dict{Int64, Vector{Int64}}()
    blockDurations = Dict{Int64, Float64}()
    delayIDs_tmp = Dict{Int64, Float64}()
    while true
        if version_combined <= 1002001
            NumberBlockEvents = 7
        else
            NumberBlockEvents = 8
        end

        read_event = readline(io)
        !isempty(read_event) || break
        blockEvents = parse.(Int64, split(read_event))

        if blockEvents[1] != 0
            if version_combined <= 1002001
                eventTable[blockEvents[1]] = Int64[0; blockEvents[3:end]...; 0]
            else
                eventTable[blockEvents[1]] = Int64[0; blockEvents[3:end]...]
            end

            if version_combined >= 1004000
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
    version_combined = 0
    version_major = 0
    version_minor = 0
    gradLibrary = Dict()
    def = Dict()
    signature = ""
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
                signature = read_signature(io)
            else
            end

        end
    end
    # fix blocks, gradients and RF objects imported from older versions
    if version_combined < 1004000
        # scan through the RF objects
        for i = 0:length(rfLibrary)-1
            rfLibrary[i]["data"] = [rfLibrary[i]["data"][1:3]' 0.0 rfLibrary[i]["data"][4:end]']
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
                gradLibrary[i]["data"] = [gradLibrary[i]["data"][1:2]; 0.0; gradLibrary[i]["data"][3:end]]
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
        seq += get_block(obj, i)
    end
    # Add first and last points for gradients #320
    fix_first_last_grads!(seq)
    # Final details
    # Hack for including extension and triggers, this will be done properly for #308 and #323
    seq.DEF["additional_text"] = read_Extension(extensionLibrary, triggerLibrary) #Temporary hack
    seq.DEF = merge(obj["definitions"], seq.DEF)
    # Koma specific details for reconstrucion
    seq.DEF["FileName"] = basename(filename)
    seq.DEF["PulseqVersion"] = version_combined
    seq.DEF["signature"] = signature
    # Guessing recon dimensions
    # ideally all estimates of recon dimensions in one place, as late as possible, i.e. not here, non-Cartesean ?? *** CAC 240708
    #seq.DEF["Ns"] = get(seq.DEF, "Ns", length(unique(seq.RF.Δf))) #slices or slabs
    #seq.DEF["Nx"] = get(seq.DEF, "Nx", maximum(adc.N for adc = seq.ADC)) #readout length
    #seq.DEF["Ny"] = get(seq.DEF, "Ny", sum(map(is_ADC_on, seq)) ÷ seq.DEF["Nx"]) #pe1
    #seq.DEF["Nz"] = get(seq.DEF, "Nz", sum(map(is_ADC_on, seq)) ÷ seq.DEF["Ny"]) #pe2
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
    G = Grad(0.0,0.0)
    if isempty(gradLibrary) || i==0
        return G
    end

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
            gT = diff(gt) * Δt_gr
            G = Grad(gA,gT,0.0,0.0,delay)
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

    if isempty(rfLibrary) || i==0
        return reshape([RF(0.0,0.0)], 1, 1)
    end

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
        rfA = decompress_shape(shapeLibrary[mag_id]...)
        rfϕ = decompress_shape(shapeLibrary[phase_id]...)
        @assert all(rfϕ.>=0) "[RF id $i] Phase waveform rfϕ must have non-negative samples (1.>=rfϕ.>=0). "
        Nrf = shapeLibrary[mag_id][1] - 1
        rfAϕ = amplitude .* rfA .* exp.(1im*(2π*rfϕ .+ phase))
    else
        rfAϕ = ComplexF64[0.0]
    end
    #Creating timings
    if time_shape_id == 0 #no time waveform
        rfT = Nrf * Δt_rf
    else
        rft = decompress_shape(shapeLibrary[time_shape_id]...)
        rfT = diff(rft) * Δt_rf
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

    if isempty(adcLibrary) || i==0
        return [ADC(0, 0)]
    end

    #Unpacking
    #(1)num (2)dwell (3)delay (4)freq (5)phase
    a = adcLibrary[i]["data"]
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
    D = Float64[max(obj["blockDurations"][i], dur(Gx), dur(Gy), dur(Gz), dur(R[1]), dur(A[1]))]

    #Extensions
    E = Dict{String, Any}()#read_Extension(obj["extensionLibrary"], iext, i)

    #Sequence block definition
    s = Sequence(G,R,A,D,E)
    s
end

#######################################################################
#######################################################################
#######################################################################
# @EventLibrary/EventLibrary.m

mutable struct EventLibrary
    keys::Vector{Int}
    data::Vector{Dict{Symbol, Any}}
    lengths::Vector{Int}
    type::Vector{Char}
    keymap::Dict{String, Int}
    next_free_id::Int
    # id_hit_count::Vector{Int}  # Uncomment if needed

    function EventLibrary()
        new(keys=[], data=Dict{Symbol, Any}[], lengths=[], type=Char[], keymap=Dict{String, Int}(), next_free_id=1)
        # id_hit_count=[]  # Uncomment if needed
    end
end

function find(obj::EventLibrary, data)
    data_string = join(data, " ")
    try
        id = obj.keymap[data_string]
        found = true
        # obj.id_hit_count[id] += 1  # Uncomment if needed
    catch
        id = obj.next_free_id
        found = false
    end
    return id, found
end

function find_or_insert(obj::EventLibrary, data, type)
    data_string = join(data, " ")
    try
        id = obj.keymap[data_string]
        found = true
        # obj.id_hit_count[id] += 1  # Uncomment if needed
    catch
        id = obj.next_free_id
        found = false
        obj.keys = vcat(obj.keys, [id])
        obj.data = vcat(obj.data, Dict(:array => data))
        obj.lengths = vcat(obj.lengths, [length(data)])
        if length(type) > 0
            obj.type = vcat(obj.type, type)
        end
        obj.keymap[data_string] = id
        obj.next_free_id = id + 1
    end
    return id, found
end

function insert(obj::EventLibrary, id, data, type=Char[])
    if id == 0
        id = obj.next_free_id
    end
    obj.keys = vcat(obj.keys, [id])
    obj.data = vcat(obj.data, Dict(:array => data))
    obj.lengths = vcat(obj.lengths, [length(data)])
    if length(type) > 0
        obj.type = vcat(obj.type, type)
    end
    data_string = join(data, " ")
    obj.keymap[data_string] = id
    if id >= obj.next_free_id
        obj.next_free_id = id + 1
    end
end

function update(obj::EventLibrary, id, old_data, new_data, type=Char[])
    if length(obj.keys) >= id
        data_string = join(old_data, " ")
        delete!(obj.keymap, data_string)
    end
    if length(type) > 0
        insert(obj, id, new_data, type)
    else
        insert(obj, id, new_data)
    end
end

function update_data(obj::EventLibrary, id, old_data, new_data, type=Char[])
    # (id, found) = find(obj, old_data)
    # if found
    if length(type) > 0
        update(obj, id, old_data, new_data, type)
    else
        update(obj, id, old_data, new_data)
    end
    # else
    #     if length(type) > 0
    #         insert(obj, id, new_data, type)
    #     else
    #         insert(obj, id, new_data)
    #     end
    # end
end

function get(obj::EventLibrary, id)
    out = Dict(:key => obj.keys[id], :data => obj.data[id][:array], :length => obj.lengths[id], :type => obj.type[id])
    return out
end

#######################################################################
#######################################################################
#######################################################################
# @EventLibrary/find_mat.m

function find_mat(keys, data, lengths, newData)
    # find_mat Lookup a data structure in the given library.
    #   idx=find_mat(keys, data, lengths, newDat) Return the index of the
    #   newDat in the library. If the data doesn't exist in the library
    #   then the index for the next new entry is returned.
    #
    #   This function is compatible with the MATLAB Coder to generate a mex
    #   file for faster execution.
    #
    #   See also EventLibrary

    found = false
    id = 0

    for i in 1:length(data)
        if lengths[i] == length(newData) && norm(data[i][:array] .- newData) < 1e-6
            id = keys[i]
            found = true
            break
        end
    end

    if isempty(keys)
        id = 1
    elseif ~found
        id = maximum(keys) + 1
    end

    return id, found
end

############################################################################################

mutable struct SEQ
    version_major::Int
    version_minor::Int
    version_revision::Int
    rfRasterTime::Float64
    gradRasterTime::Float64
    adcRasterTime::Float64
    blockDurationRaster::Float64
    definitions::Dict{String, Any}

    blockEvents::Vector{Any}
    blockDurations::Vector{Any}
    rfLibrary::EventLibrary
    gradLibrary::EventLibrary
    adcLibrary::EventLibrary
    trigLibrary::EventLibrary
    labelsetLibrary::EventLibrary
    labelincLibrary::EventLibrary
    extensionLibrary::EventLibrary
    shapeLibrary::EventLibrary
    extensionStringIDs::Vector{String}
    extensionNumericIDs::Vector{Int}

    signatureType::String
    signatureFile::String
    signatureValue::String
    sys
end

function SEQ(; sys=nothing)
    version_major = 1
    version_minor = 4
    version_revision = 1
    definitions = Dict{String, Any}()
    gradLibrary = EventLibrary()
    shapeLibrary = EventLibrary()
    rfLibrary = EventLibrary()
    adcLibrary = EventLibrary()
    trigLibrary = EventLibrary()
    labelsetLibrary = EventLibrary()
    labelincLibrary = EventLibrary()
    extensionLibrary = EventLibrary()
    extensionStringIDs = []
    extensionNumericIDs = []
    blockEvents = []

    if isnothing(sys)
        sys = opts()
    end

    rfRasterTime = sys.rfRasterTime
    gradRasterTime = sys.gradRasterTime
    adcRasterTime = sys.adcRasterTime
    blockDurationRaster = sys.blockDurationRaster

    setDefinition(key, val) = (definitions[key] = val)
    getDefinition(key) = get(definitions, key, nothing)

    signatureType = ""
    signatureFile = ""
    signatureValue = ""

    return SEQ(
        version_major,
        version_minor,
        version_revision,
        rfRasterTime,
        gradRasterTime,
        adcRasterTime,
        blockDurationRaster,
        definitions,
        blockEvents,
        blockDurations,
        rfLibrary,
        gradLibrary,
        adcLibrary,
        trigLibrary,
        labelsetLibrary,
        labelincLibrary,
        extensionLibrary,
        shapeLibrary,
        extensionStringIDs,
        extensionNumericIDs,
        signatureType,
        signatureFile,
        signatureValue,
        sys
    )
end

###############################################
function duration(obj::SEQ)
    """
    duration()
        Returns the total duration of the sequence
        optionally returns the total count of events
    """

    # Loop over blocks and gather statistics
    numBlocks = length(obj.blockEvents)
    eventCount = (numBlocks > 0 && nargout() > 2) ? zeros(size(obj.blockEvents[1])) : nothing
    duration = 0

    for iB in 1:numBlocks
        if nargout() > 2
            eventCount .= eventCount .+ (obj.blockEvents[iB] .> 0)
        end
        duration += obj.blockDurations[iB]
    end

    return duration, numBlocks, eventCount
end


function checkTiming(obj::SEQ)
    """
    checkTiming()
        Checks timing of all blocks and objects in the sequence
        optionally returns the detailed error log as a cell array
        of strings. This function also modifies the sequence
        object by adding the field "TotalDuration" to sequence
        definitions.
    """

    # Loop over blocks and gather statistics
    numBlocks = length(obj.blockEvents)
    is_ok = true
    errorReport = []
    totalDuration = 0

    for iB in 1:numBlocks
        b = obj.getBlock(iB)
        ind = ~structfun(@isempty, b)
        fn = fieldnames(b)
        ev = [b.(f) for f in fn[ind]]

        res, rep, dur = mr.checkTiming(obj.sys, ev...)

        is_ok = is_ok && res

        # check the stored block duration
        if abs(dur - obj.blockDurations[iB]) > eps
            rep *= "inconsistency between the stored block duration and the duration of the block content"
            is_ok = false
            dur = obj.blockDurations[iB]
        end

        # check that block duration is aligned to the blockDurationRaster
        bd = obj.blockDurations[iB] / obj.blockDurationRaster
        bdr = round(bd)
        if abs(bdr - bd) >= 1e-6
            rep *= "block duration is not aligned to the blockDurationRaster"
            is_ok = false
        end

        # check RF dead times
        if ~isempty(b.rf)
            if b.rf.delay - b.rf.deadTime < -eps
                rep *= "delay of $(b.rf.delay*1e6)us is smaller than the RF dead time $(b.rf.deadTime*1e6)us"
                is_ok = false
            end
            if b.rf.delay + b.rf.t[end] + b.rf.ringdownTime - dur > eps
                rep *= "time between the end of the RF pulse at $(b.rf.delay+b.rf.t[end])*1e6 and the end of the block at $(dur*1e6)us is shorter than rfRingdownTime"
                is_ok = false
            end
        end

        # check ADC dead times
        if ~isempty(b.adc)
            if b.adc.delay - obj.sys.adcDeadTime < -eps
                rep *= " adc.delay < system.adcDeadTime"
                is_ok = false
            end
            if b.adc.delay + b.adc.numSamples * b.adc.dwell + obj.sys.adcDeadTime - dur > eps
                rep *= " adc: system.adcDeadTime (post-adc) violation"
                is_ok = false
            end
        end

        # update report
        if ~isempty(rep)
            push!(errorReport, "   Block: $iB $rep\n")
        end

        totalDuration += dur
    end

    # check whether all gradients in the last block are ramped down properly
    if ~isempty(ev) && isstruct(ev)
        for en in 1:length(ev)
            if length(ev[en]) == 1 && ev[en].type == "grad"  # length(ev{en})==1 excludes arrays of extensions
                if ev[en].last != 0  # must be > sys.slewRate*sys.gradRasterTime
                    push!(errorReport, "   Block: $iB gradients do not ramp to 0 at the end of the sequence\n")
                end
            end
        end
    end

    obj.setDefinition("TotalDuration", totalDuration)
    return is_ok, errorReport
end

function getDefinition(obj::SEQ, key)
    """
    getDefinition Return the values of custom definition.
    val = getDefinitions(seqObj, key) Return value of the
    definition specified by the key.

    These definitions can be added manually or read from the
    header of a sequence file defined in the sequence header.
    An empty array is returned if the key is not defined.

    See also setDefinition
    """
    value = haskey(obj.definitions, key) ? obj.definitions[key] : []
    return value
end

function setDefinition(obj::SEQ, key, val)
    """
    setDefinition Modify a custom definition of the sequence.
    setDefinition(obj, def, val) Set the user definition 'key'
    to value 'val'. If the definition does not exist, it will be
    created.

    See also getDefinition
    """
    if key == "FOV"
        # Issue a warning if FOV is too large, e.g., is in mm
        if maximum(val) > 1
            println("WARNING: definition FOV uses values exceeding 1m. New Pulseq interpreters expect values in units of meters!")
        end
    end
    obj.definitions[key] = val
end

function addBlock(obj::SEQ, varargin...)
    """
    addBlock Add a new block to the sequence.
    addBlock(obj, blockStruct) Adds a sequence block with
    provided as a block structure

    addBlock(obj, e1, e2, ...) Adds a block with multiple
    events e1, e2, etc.

    See also  setBlock, makeAdc, makeTrapezoid, makeSincPulse
    """
    setBlock(obj, length(obj.blockEvents) + 1, varargin...)
end

function modGradAxis(obj::SEQ, axis, modifier)
    # modGradAxis Invert or scale all gradients along the corresponding axis/channel.
    # The function acts on all gradient objects already added to the sequence object

    channelNum = findfirst(x -> x in ["x", "y", "z"], axis)
    otherChans = findall(x -> !(x in ["x", "y", "z"]), axis)

    # go through all event table entries and list gradient objects in the library
    paren2 = (x, varargin) -> x[:, varargin]
    allGradEvents = paren2(vertcat(obj.blockEvents...), 3:5)

    selectedEvents = Set(allGradEvents[:, channelNum])
    selectedEvents = Set(filter(x -> x != 0, selectedEvents))  # eliminate 0

    otherEvents = Set(allGradEvents[:, otherChans])
    @assert isempty(intersect(selectedEvents, otherEvents),
                   "ERROR: the same gradient event is used on multiple axes, this is not yet supported by modGradAxis()")

    for i in selectedEvents
        # based on the above, we just patch the first element of the gradient library data entries
        obj.gradLibrary.data[i].array[1] *= modifier

        if obj.gradLibrary.type[i] == "g" && obj.gradLibrary.lengths[i] == 5
            # need to update .first and .last fields
            obj.gradLibrary.data[i].array[4] *= modifier
            obj.gradLibrary.data[i].array[5] *= modifier
        end
    end
end

function flipGradAxis(obj::SEQ, axis)
    # flipGradAxis Invert all gradients along the corresponding axis/channel.
    # The function acts on all gradient objects already added to the sequence object
    modGradAxis(obj, axis, -1)
end

function rf = rfFromLibData(obj::SEQ, libData, use)
    rf = Dict()

    rf["type"] = "rf"

    amplitude = libData[1]
    magShape = libData[2]
    phaseShape = libData[3]
    shapeData = obj.shapeLibrary.data[magShape].array
    compressed.num_samples = shapeData[1]
    compressed.data = shapeData[2:end]
    mag = mr.decompressShape(compressed)
    shapeData = obj.shapeLibrary.data[phaseShape].array
    compressed.num_samples = shapeData[1]
    compressed.data = shapeData[2:end]
    phase = mr.decompressShape(compressed)
    rf["signal"] = amplitude .* mag .* exp(1j * 2 * pi * phase)
    timeShape = libData[4]
    if timeShape > 0
        shapeData = obj.shapeLibrary.data[timeShape].array
        compressed.num_samples = shapeData[1]
        compressed.data = shapeData[2:end]
        rf["t"] = mr.decompressShape(compressed) .* obj.rfRasterTime
        rf["shape_dur"] = ceil((rf["t"][end] - eps()) / obj.rfRasterTime) * obj.rfRasterTime
    else
        # generate default time raster on the fly
        rf["t"] = ((1:length(rf["signal"])) .- 0.5) .* obj.rfRasterTime
        rf["shape_dur"] = length(rf["signal"]) * obj.rfRasterTime
    end

    rf["delay"] = libData[5]
    rf["freqOffset"] = libData[6]
    rf["phaseOffset"] = libData[7]

    rf["deadTime"] = obj.sys.rfDeadTime
    rf["ringdownTime"] = obj.sys.rfRingdownTime

    if length(libData) < 8
        libData[8] = 0
    end
    rf["deadTime"] = libData[9]

    if length(libData) < 9
        libData[9] = 0
    end
    rf["ringdownTime"] = libData[9]

    if nargin > 2
        rf["use"] = begin
            if use == "e"
                return "excitation"
            elseif use == "r"
                return "refocusing"
            elseif use == "i"
                return "inversion"
            elseif use == "s"
                return "saturation"
            elseif use == "p"
                return "preparation"
            else
                return "undefined"
            end
        end
    else
        rf["use"] = "undefined"
    end
end


function registerRfEvent(obj::SEQ, event)
    # registerRfEvent Add the event to the libraries (object,
    # shapes, etc and return the event's ID. This ID can be stored in
    # the object to accelerate addBlock()

    mag = abs(event["signal"])
    amplitude = maximum(mag)
    mag = mag / amplitude
    phase = angle.(event["signal"])
    phase[phase .< 0] .+= 2 * π
    phase = phase / (2 * π)
    may_exist = true

    if haskey(event, "shapeIDs")
        shapeIDs = event["shapeIDs"]
    else
        shapeIDs = [0, 0, 0]

        magShape = mr.compressShape(mag[:])
        data = [magShape.num_samples; magShape.data]
        shapeIDs[1], found = obj.shapeLibrary.find_or_insert(data)
        may_exist = may_exist && found

        phaseShape = mr.compressShape(phase)
        data = [phaseShape.num_samples; phaseShape.data]
        shapeIDs[2], found = obj.shapeLibrary.find_or_insert(data)
        may_exist = may_exist && found

        timeShape = mr.compressShape(event["t"] / obj.rfRasterTime)
        # time shape is stored in units of RF raster
        if length(timeShape.data) == 4 && all(timeShape.data .== [0.5, 1, 1, timeShape.num_samples - 3])
            shapeIDs[3] = 0
        else
            data = [timeShape.num_samples; timeShape.data]
            shapeIDs[3], found = obj.shapeLibrary.find_or_insert(data)
            may_exist = may_exist && found
        end
    end

    use = "u"
    if haskey(event, "use")
        use = begin
            if event["use"] == "excitation"
                return "e"
            elseif event["use"] == "refocusing"
                return "r"
            elseif event["use"] == "saturation"
                return "s"
            elseif event["use"] == "preparation"
                return "p"
            else
                return "u"
            end
        end
    end

    data = [amplitude; shapeIDs[1]; shapeIDs[2]; shapeIDs[3];
            event["delay"]; event["freqOffset"]; event["phaseOffset"]]  # ;
            #event["deadTime"]; event["ringdownTime"]]

    if may_exist
        id, _ = obj.rfLibrary.find_or_insert(data, use)
    else
        id = obj.rfLibrary.insert(0, data, use)
    end

    return id, shapeIDs
end

function registerGradEvent(obj::SEQ, event)
    may_exist = true

    id = 0
    shapeIDs = [0, 0]

    if event.type == "grad"
        amplitude = maximum(abs.(event.waveform))
        if amplitude > 0
            _, _, fnz = findfirst(x -> x != 0, event.waveform)  # find the first non-zero value and make it positive
            amplitude *= sign(fnz)
        end

        if hasfield(event, :shapeIDs)
            shapeIDs = event.shapeIDs
        else
            shapeIDs = [0, 0]

            if amplitude != 0
                g = event.waveform / amplitude
            else
                g = event.waveform
            end

            c_shape = compressShape(g)
            s_data = [c_shape.num_samples; c_shape.data]
            shapeIDs[1], found = find_or_insert(obj.shapeLibrary, s_data...)
            may_exist &= found

            c_time = compressShape(event.tt / obj.gradRasterTime)
            if !(length(c_time.data) == 4 && all(c_time.data .== [0.5, 1, 1, c_time.num_samples - 3]))
                t_data = [c_time.num_samples; c_time.data]
                shapeIDs[2], found = find_or_insert(obj.shapeLibrary, t_data...)
                may_exist &= found
            end
        end

        data = [amplitude; shapeIDs; event.delay; event.first; event.last]
    elseif event.type == "trap"
        data = [event.amplitude; event.riseTime; event.flatTime; event.fallTime; event.delay]
    else
        error("unknown gradient type passed to registerGradEvent()")
    end

    if may_exist
        id, found = find_or_insert(obj.gradLibrary, data, event.type[1])
        if !found
            id = 0
        end
    else
        id = insert(obj.gradLibrary, 0, data, event.type[1])
    end

    return id, shapeIDs
end

function registerAdcEvent(obj::SEQ, event)
    data = [event.numSamples, event.dwell, max(event.delay, event.deadTime),
            event.freqOffset, event.phaseOffset, event.deadTime]

    id, found = find_or_insert(obj.adcLibrary, data...)

    if !found
        id = 0
    end

    return id
end

function registerControlEvent(obj::SEQ, event)
    event_type = findfirst(x -> x == event.type, ["output", "trigger"])

    if event_type == 1
        event_channel = findfirst(x -> x == event.channel, ["osc0", "osc1", "ext1"])
    elseif event_type == 2
        event_channel = findfirst(x -> x == event.channel, ["physio1", "physio2"])
    else
        error("unsupported control event type")
    end

    data = [event_type, event_channel, event.delay, event.duration]
    id, found = find_or_insert(obj.trigLibrary, data...)

    if !found
        id = 0
    end

    return id
end

function registerLabelEvent(obj::SEQ, event)
    label_id = findfirst(x -> x == event.label, mr.getSupportedLabels())
    data = [event.value, label_id]

    id, found = 0, false

    if event.type == "labelset"
        id, found = find_or_insert(obj.labelsetLibrary, data...)
    elseif event.type == "labelinc"
        id, found = find_or_insert(obj.labelincLibrary, data...)
    else
        error("unknown label type passed to registerLabelEvent()")
    end

    if !found
        id = 0
    end

    return id
end

# TODO: Replacing blocks in the middle of sequence can cause unused
# events in the libraries. These can be detected and pruned.
function setBlock(obj, index, varargin)
    # setBlock Replace or add sequence block.
    # setBlock(obj, index, bStruct) Replace block at index with new
    # block provided as block structure.
    #
    # setBlock(obj, index, e1, e2, ...) Create a new block from
    # events and store at position given by index.
    #
    # The block or events are provided in uncompressed form and
    # will be stored in the compressed, non-redundant internal
    # libraries.
    #
    # See also  getBlock, addBlock

    # Convert block structure to cell array of events
    varargin = mr.block2events(varargin)

    obj.blockEvents[index] .= zeros(7)
    duration = 0

    check_g = Dict{Any, Any}()  # Dictionary containing a structure, each with the index and pairs of gradients/times
    extensions = []

    # Loop over events adding to library if necessary and creating block event structure.
    for event in varargin
        event_type = event.type

        if event_type == "rf"
            id = haskey(event, "id") ? event.id : obj.registerRfEvent(event)
            obj.blockEvents[index][2] = id
            duration = max(duration, event.shape_dur + event.delay + event.ringdownTime)

        elseif event_type == "grad"
            channelNum = findfirst(x -> x == event.channel, ["x", "y", "z"])
            idx = 2 + channelNum

            grad_start = event.delay + floor(event.tt[1] / obj.gradRasterTime + 1e-10) * obj.gradRasterTime
            grad_duration = event.delay + ceil(event.tt[end] / obj.gradRasterTime - 1e-10) * obj.gradRasterTime

            check_g[channelNum] = Dict("idx" => idx, "start" => [grad_start, event.first], "stop" => [grad_duration, event.last])

            id = haskey(event, "id") ? event.id : obj.registerGradEvent(event)
            obj.blockEvents[index][idx] = id
            duration = max(duration, grad_duration)

        elseif event_type == "trap"
            channelNum = findfirst(x -> x == event.channel, ["x", "y", "z"])
            idx = 2 + channelNum

            check_g[channelNum] = Dict("idx" => idx, "start" => [0, 0], "stop" => [event.delay + event.riseTime + event.fallTime + event.flatTime, 0])

            id = haskey(event, "id") ? event.id : obj.registerGradEvent(event)
            obj.blockEvents[index][idx] = id
            duration = max(duration, event.delay + event.riseTime + event.flatTime + event.fallTime)

        elseif event_type == "adc"
            id = haskey(event, "id") ? event.id : obj.registerAdcEvent(event)
            obj.blockEvents[index][6] = id
            duration = max(duration, event.delay + event.numSamples * event.dwell + event.deadTime)

        elseif event_type == "delay"
            duration = max(duration, event.delay)

        elseif event_type in ["output", "trigger"]
            id = haskey(event, "id") ? event.id : obj.registerControlEvent(event)
            ext = Dict("type" => obj.getExtensionTypeID("TRIGGERS"), "ref" => id)
            push!(extensions, ext)
            duration = max(duration, event.delay + event.duration)

        elseif event_type in ["labelset", "labelinc"]
            id = haskey(event, "id") ? event.id : obj.registerLabelEvent(event)
            ext = Dict("type" => obj.getExtensionTypeID(uppercase(event.type)), "ref" => id)
            push!(extensions, ext)
        end
    end

    if !isempty(extensions)
        sort!(extensions, by = x -> x["ref"])
        all_found = true
        id = 0

        for ext in extensions
            data = [ext["type"], ext["ref"], id]
            id, found = obj.extensionLibrary.find(data)
            all_found = all_found && found

            if !found
                break
            end
        end

        if !all_found
            for ext in extensions
                data = [ext["type"], ext["ref"], id]
                id, found = obj.extensionLibrary.find(data)

                if !found
                    obj.extensionLibrary.insert(id, data)
                end
            end
        end

        obj.blockEvents[index][7] = id
    end

    # Check if connection to the previous block is correct using check_g
    for (channelNum, cg) in check_g
        if isempty(cg)
            continue
        end

        # Check the start
        if abs(cg["start"][2]) > obj.sys.maxSlew * obj.sys.gradRasterTime
            if cg["start"][1] != 0
                error("Error in block $index: No delay allowed for gradients which start with a non-zero amplitude")
            end

            if index > 1
                prev_nonempty_block = findlast(x -> obj.blockDurations[x] > 0, 1:length(obj.blockDurations))
                prev_id = obj.blockEvents[prev_nonempty_block][cg["idx"]]

                if prev_id != 0
                    prev_lib = obj.gradLibrary.get(prev_id)

                    prev_dat = prev_lib.data
                    prev_type = prev_lib.type

                    if prev_type == "t"
                        error("Error in block $index: Two consecutive gradients need to have the same amplitude at the connection point, this is not possible if the previous gradient is a simple trapezoid")
                    elseif prev_type == "g"
                        last = prev_dat[6]
                        if abs(last - cg["start"][2]) > obj.sys.maxSlew * obj.sys.gradRasterTime
                            error("Error in block $index: Two consecutive gradients need to have the same amplitude at the connection point")
                        end
                    end
                else
                    error("Error in block $index: Gradient starting at non-zero value need to be preceded by a compatible gradient")
                end
            else
                error("First gradient in the first block has to start at 0.")
            end
        end

        # Check if gradients, which do not end at 0, are as long as the block itself.
        if cg["stop"][2] > obj.sys.maxSlew * obj.sys.gradRasterTime && abs(cg["stop"][1] - duration) > 1e-7
            error("Error in block $index: A gradient that doesn't end at zero needs to be aligned to the block boundary")
        end
    end

    obj.blockDurations[index] = duration
end




############################################################################################

function getBlock(obj, index)
    block = (
        blockDuration = 0,
        rf = [],
        gx = [],
        gy = [],
        gz = [],
        adc = [],
        trig = [],
        label = []
    )

    block.rf = []
    eventInd = obj.blockEvents[index]

    if eventInd[7] > 0
        # we have extensions -- triggers, labels, etc
        # we will eventually isolate this into a separate function
        nextExtID = eventInd[7]
        while nextExtID != 0
            extData = obj.extensionLibrary.data[nextExtID].array
            # format: extType, extID, nextExtID
            extTypeStr = getExtensionTypeString(obj, extData[1])
            if extTypeStr == "TRIGGERS"
                trigger_types = ["output", "trigger"]
                data = obj.trigLibrary.data[extData[2]].array
                trig.type = trigger_types[data[1]]
                if data[1] == 1
                    trigger_channels = ["osc0", "osc1", "ext1"]
                    trig.channel = trigger_channels[data[2]]
                elseif data[1] == 2
                    trigger_channels = ["physio1", "physio2"]
                    trig.channel = trigger_channels[data[2]]
                else
                    error("unsupported trigger event type")
                end
                trig.delay = data[3]
                trig.duration = data[4]
                # allow for multiple triggers per block
                if hasfield(block, "trig")
                    push!(block.trig, trig)
                else
                    block.trig = [trig]
                end
            elseif extTypeStr == ("LABELSET", "LABELINC")
                label.type = lowercase(extTypeStr)
                supported_labels = mr.getSupportedLabels()
                if extTypeStr == "LABELSET"
                    data = obj.labelsetLibrary.data[extData[2]].array
                else
                    data = obj.labelincLibrary.data[extData[2]].array
                end
                label.label = supported_labels[data[2]]
                label.value = data[1]
                # allow for multiple labels per block
                if hasfield(block, "label")
                    push!(block.label, label)
                else
                    block.label = [label]
                end
            else
                error("unknown extension ID $(extData[1])")
            end
            # now update nextExtID
            nextExtID = extData[3]
        end
    end

    if eventInd[2] > 0
        if length(obj.rfLibrary.type) >= eventInd[2]
            block.rf = obj.rfFromLibData(obj.rfLibrary.data[eventInd[2]].array, obj.rfLibrary.type[eventInd[2]])
        else
            block.rf = obj.rfFromLibData(obj.rfLibrary.data[eventInd[2]].array)  # undefined type/use
        end
    end

    gradChannels = ["gx", "gy", "gz"]
    for (i, gradChannel) in enumerate(gradChannels)
        if eventInd[2 + i] > 0
            type = obj.gradLibrary.type[eventInd[2 + i]]
            libData = obj.gradLibrary.data[eventInd[2 + i]].array
            grad.type = ifelse(type == "t", "trap", "grad")
            grad.channel = gradChannel[2]
            if grad.type == "grad"
                amplitude = libData[1]
                shapeId = libData[2]
                timeId = libData[3]
                delay = libData[4]
                shapeData = obj.shapeLibrary.data[shapeId].array
                compressed.num_samples = shapeData[1]
                compressed.data = shapeData[2:end]
                try
                    g = mr.decompressShape(compressed)
                catch
                    @printf("  mr.decompressShape() failed for shapeId %d\n", shapeId)
                    error("mr.decompressShape() failed for shapeId %d", shapeId)
                end
                grad.waveform = amplitude * g
                if timeId == 0
                    grad.tt = ((1:length(g)) .- 0.5) * obj.gradRasterTime  # TODO: eventually we may remove these true-times
                    t_end = length(g) * obj.gradRasterTime
                else
                    tShapeData = obj.shapeLibrary.data[timeId].array
                    compressed.num_samples = tShapeData[1]
                    compressed.data = tShapeData[2:end]
                    try
                        grad.tt = mr.decompressShape(compressed) * obj.gradRasterTime
                    catch
                        @printf("  mr.decompressShape() failed for shapeId %d\n", shapeId)
                        error("mr.decompressShape() failed for shapeId %d", shapeId)
                    end
                    @assert length(grad.waveform) == length(grad.tt)
                    t_end = grad.tt[end]
                end
                grad.shape_id = shapeId  # needed for the second pass of read()
                grad.time_id = timeId  # needed for the second pass of read()
                grad.delay = delay
                grad.shape_dur = t_end
                if length(libData) > 5
                    grad.first = libData[5]
                    grad.last = libData[6]
                else
                    if hasfield(grad, "first"), grad = rmfield(grad, "first"); end
                    if hasfield(grad, "last"), grad = rmfield(grad, "last"); end
                end
            else
                grad.amplitude = libData[1]
                grad.riseTime = libData[2]
                grad.flatTime = libData[3]
                grad.fallTime = libData[4]
                grad.delay = libData[5]
                grad.area = grad.amplitude * (grad.flatTime + grad.riseTime / 2 + grad.fallTime / 2)
                grad.flatArea = grad.amplitude * grad.flatTime
            end

            block.(gradChannel) = grad
        end
    end

    if eventInd[6] > 0
        libData = obj.adcLibrary.data[eventInd[6]].array
        if length(libData) < 6
            libData = vcat(libData, 0)
        end
        adc = cell2struct(num2cell(libData),
                          [:numSamples, :dwell, :delay, :freqOffset, :phaseOffset, :deadTime], 2)
        adc.type = "adc"
        block.adc = adc
    end

    block.blockDuration = obj.blockDurations[index]
    return block
end

function calculateKspaceUnfunc(obj; trajectory_delay=0)
    # calculate the k-space trajectory of the entire pulse sequence
    # optional parameter 'trajectory_delay' sets the compensation factor to align ADC and gradients in the reconstruction
    # Return values: ktraj_adc, ktraj, t_excitation, t_refocusing

    if abs(trajectory_delay) > 100e-6
        println("Warning: trajectory delay of $(trajectory_delay * 1e6) us is suspiciously high")
    end

    # initialise the counters and accumulator objects
    c_excitation = 0
    c_refocusing = 0
    c_adcSamples = 0

    # loop through the blocks to prepare preallocations
    for iB = 1:length(obj.blockEvents)
        block = getBlock(obj, iB)
        if !isempty(block.rf)
            if !hasfield(block.rf, "use") || block.rf.use == "excitation" || block.rf.use == "undefined"
                c_excitation += 1
            elseif block.rf.use == "refocusing"
                c_refocusing += 1
            end
        end
        if !isempty(block.adc)
            c_adcSamples += block.adc.numSamples
        end
    end

    # preallocate arrays
    t_excitation = zeros(c_excitation)
    t_refocusing = zeros(c_refocusing)
    ktime = zeros(c_adcSamples)
    current_dur = 0
    c_excitation = 1
    c_refocusing = 1
    kcouter = 1
    traj_recon_delay = trajectory_delay

    # go through the blocks and collect RF and ADC timing data
    for iB = 1:length(obj.blockEvents)
        block = getBlock(obj, iB)
        if !isempty(block.rf)
            rf = block.rf
            t = rf.delay + calcRfCenter(rf)
            if !hasfield(block.rf, "use") || block.rf.use == "excitation" || block.rf.use == "undefined"
                t_excitation[c_excitation] = current_dur + t
                c_excitation += 1
            elseif block.rf.use == "refocusing"
                t_refocusing[c_refocusing] = current_dur + t
                c_refocusing += 1
            end
        end
        if !isempty(block.adc)
            ktime[kcouter:kcouter + block.adc.numSamples - 1] = ((0:block.adc.numSamples - 1) .+ 0.5) * block.adc.dwell + block.adc.delay + current_dur + traj_recon_delay
            kcouter += block.adc.numSamples
        end
        current_dur += obj.blockDurations[iB]
    end

    # calculate the actual k-space trajectory based on the gradient waveforms
    gw = gradient_waveforms(obj)
    i_excitation = round.(Int, t_excitation / obj.gradRasterTime)
    i_refocusing = round.(Int, t_refocusing / obj.gradRasterTime)

    i_periods = sort([1; i_excitation .+ 1; i_refocusing .+ 1; size(gw, 2) + 1])
    ii_next_excitation = min(length(i_excitation), 1)
    ii_next_refocusing = min(length(i_refocusing), 1)
    ktraj = zeros(size(gw))
    k = [0.0, 0.0, 0.0]

    for i in 1:(length(i_periods) - 1)
        i_period_end = i_periods[i + 1] - 1
        k_period = cumsum([k gw[:, i_periods[i]:i_period_end] * obj.gradRasterTime], dims=2)
        ktraj[:, i_periods[i]:i_period_end] = k_period[:, 2:end]
        k = k_period[:, end]

        if ii_next_excitation > 0 && i_excitation[ii_next_excitation] == i_period_end
            fill!(k, 0.0)
            ktraj[:, i_period_end] .= NaN  # we use NaNs to mark the excitation point, they interrupt the plots
            ii_next_excitation = min(length(i_excitation), ii_next_excitation + 1)
        end

        if ii_next_refocusing > 0 && i_refocusing[ii_next_refocusing] == i_period_end
            k .= -k
            ii_next_refocusing = min(length(i_refocusing), ii_next_refocusing + 1)
        end
    end

    # calculate the k-space positions at the ADC time points
    # sample the k-space positions at the ADC time points
    ktraj_adc = interpolate((1:size(ktraj, 2)) * obj.gradRasterTime, ktraj, ktime, extrapolation=Flat())[1:end-1, :]
    t_adc = ktime  # we now also return the sampling time points
    return ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc
end

############################################################################################
function gradient_waveforms1(obj)
    # gradient_waveforms()
    # Decompress the entire gradient waveform
    # Returns an array of gradient_axes x timepoints
    # gradient_axes is typically 3.

    duration, numBlocks, ~ = obj.duration()

    wave_length = ceil(Int, duration / obj.gradRasterTime)
    grad_channels = 3
    grad_waveforms = zeros(grad_channels, wave_length)
    gradChannels = ["gx", "gy", "gz"]

    t0 = 0
    t0_n = 0
    for iB = 1:numBlocks
        block = getBlock(obj, iB)
        for j = 1:length(gradChannels)
            grad = getfield(block, gradChannels[j], nothing)
            if !isempty(grad)
                if grad.type == "grad"
                    nt_start = round(Int, grad.delay / obj.gradRasterTime)
                    waveform = grad.waveform
                else
                    nt_start = round(Int, grad.delay / obj.gradRasterTime)
                    if abs(grad.flatTime) > eps
                        t = cumsum([0, grad.riseTime, grad.flatTime, grad.fallTime])
                        trapform = grad.amplitude * [0, 1, 1, 0]
                    else
                        t = cumsum([0, grad.riseTime, grad.fallTime])
                        trapform = grad.amplitude * [0, 1, 0]
                    end

                    tn = floor(Int, t[end] / obj.gradRasterTime)
                    t = vcat(t, t[end] + obj.gradRasterTime)
                    trapform = vcat(trapform, 0)

                    if abs(grad.amplitude) > eps
                        waveform = mr.pts2waveform(t, trapform, obj.gradRasterTime)
                    else
                        waveform = zeros(1, tn + 1)
                    end
                end
                if any(!isfinite, waveform)
                    println("Warning: Not all elements of the generated waveform are finite!")
                end
                grad_waveforms[j, (t0_n + 1 + nt_start):(t0_n + nt_start + length(waveform))] = waveform
            end
        end

        t0 += obj.blockDurations[iB]  # mr.calcDuration(block)
        t0_n = round(Int, t0 / obj.gradRasterTime)
    end
end


function waveforms_and_times(obj, appendRF=false, blockRange=[1, length(obj.blockEvents)])
    grad_channels = 3
    gradChannels = ["gx", "gy", "gz"]
    t0 = 0
    t0_n = 0
    numBlocks = blockRange[2] - blockRange[1] + 1

    shape_channels = appendRF ? length(gradChannels) + 1 : length(gradChannels)
    shape_pieces = Vector{Vector{Any}}(undef, shape_channels, numBlocks)

    tfp_excitation = []
    tfp_refocusing = []
    t_adc = []
    fp_adc = []
    curr_dur = 0
    iP = 0
    out_len = zeros(Int, shape_channels)

    for iBc in blockRange[1]:blockRange[2]
        block = getBlock(obj, iBc)
        iP += 1

        for j = 1:length(gradChannels)
            grad = getfield(block, gradChannels[j])

            if isempty(grad)
                continue
            end

            if grad.type == "grad"
                tt_rast = grad.tt / obj.gradRasterTime .+ 0.5
                if all(abs.(tt_rast .- (1:length(tt_rast))) .< 1e-6)
                    tt_chg, waveform_chg = restoreAdditionalShapeSamples(grad.tt, grad.waveform, grad.first, grad.last, obj.gradRasterTime, iBc)
                    out_len[j] += length(tt_chg)
                    shape_pieces[j, iP] = [curr_dur + grad.delay + tt_chg, waveform_chg]
                else
                    out_len[j] += length(grad.tt)
                    shape_pieces[j, iP] = [curr_dur + grad.delay + grad.tt', grad.waveform']
                end
            else
                if abs(grad.flatTime) > eps
                    out_len[j] += 4
                    shape_pieces[j, iP] = [
                        curr_dur + grad.delay + cumsum([0, grad.riseTime, grad.flatTime, grad.fallTime]),
                        grad.amplitude * [0, 1, 1, 0]
                    ]
                else
                    if abs(grad.riseTime) > eps && abs(grad.fallTime) > eps
                        out_len[j] += 3
                        shape_pieces[j, iP] = [
                            curr_dur + grad.delay + cumsum([0, grad.riseTime, grad.fallTime]),
                            grad.amplitude * [0, 1, 0]
                        ]
                    else
                        if abs(grad.amplitude) > eps
                            println("''empty'' gradient with non-zero magnitude detected in block $iBc")
                        end
                    end
                end
            end
        end

        if !isempty(block.rf)
            rf = block.rf
            tc = calcRfCenter(rf)
            t = rf.delay + tc

            if !hasfield(block.rf, "use") || block.rf.use == "excitation" || block.rf.use == "undefined"
                push!(tfp_excitation, [curr_dur + t, rf.freqOffset, rf.phaseOffset + rf.freqOffset * tc])
            elseif block.rf.use == "refocusing"
                push!(tfp_refocusing, [curr_dur + t, rf.freqOffset, rf.phaseOffset + rf.freqOffset * tc])
            end

            if appendRF
                pre = []
                post = []
                if abs(rf.signal[1]) > 0
                    pre = [curr_dur + rf.delay + rf.t[1] - eps, 0]
                end
                if abs(rf.signal[end]) > 0
                    post = [curr_dur + rf.delay + rf.t[end] + eps, 0]
                end

                out_len[end] += length(rf.t) + size(pre, 1) + size(post, 1)
                shape_pieces[end, iP] = [pre; curr_dur + rf.delay + rf.t'; (rf.signal .* exp(1im * (rf.phaseOffset + 2 * π * rf.freqOffset * rf.t)))'; post]
            end
        end

        if !isempty(block.adc)
            ta = block.adc.dwell * ((0:(block.adc.numSamples - 1)) .+ 0.5)
            push!(t_adc, ta + block.adc.delay + curr_dur)
            push!(fp_adc, [block.adc.freqOffset * ones(Int, block.adc.numSamples); block.adc.phaseOffset .+ block.adc.freqOffset * ta])
        end

        curr_dur += obj.blockDurations[iBc]
    end

    wave_data = Vector{Vector{Any}}(undef, shape_channels)
    for j = 1:shape_channels
        wave_data[j] = zeros(2, out_len[j])
    end

    wave_cnt = zeros(Int, shape_channels)
    curr_dur = 0

    for iP in 1:numBlocks
        for j in 1:shape_channels
            if !isempty(shape_pieces[j, iP])
                wave_data_local = shape_pieces[j, iP]
                len = size(wave_data_local, 2)

                if wave_cnt[j] != 0 && wave_data[j][1, wave_cnt[j]] + obj.gradRasterTime < wave_data_local[1, 1]
                    if wave_data[j][2, wave_cnt[j]] != 0
                        println("waveforms_and_times(): forcing ramp-down from a non-zero gradient sample on axis $j at t=$(round(1e6 * wave_data[j][1, wave_cnt[j]])) us \ncheck your sequence, some calculations are possibly wrong. If using mr.makeArbitraryGrad() consider using explicit values for 'first' and 'last' and setting them correctly.")
                        wave_data[j][:, wave_cnt[j] + 1] = [wave_data[j][1, wave_cnt[j]] + obj.gradRasterTime / 2, 0]  # this is likely to cause memory reallocations
                        wave_cnt[j] += 1
                    end

                    if wave_data_local[2, 1] != 0
                        println("waveforms_and_times(): forcing ramp-up to a non-zero gradient sample on axis $j at t=$(round(1e6 * wave_data_local[1, 1])) us \ncheck your sequence, some calculations are probably wrong.  If using mr.makeArbitraryGrad() consider using explicit values for 'first' and 'last' and setting them correctly.")
                        wave_data_local = [wave_data_local[1, 1] - obj.gradRasterTime / 2 0; wave_data_local]  # this is likely to cause memory reallocations also later on
                        len += 1
                    end
                end

                if wave_cnt[j] == 0 || wave_data[j][1, wave_cnt[j]] != wave_data_local[1, 1]
                    wave_data[j][:, wave_cnt[j] + (1:len)] = wave_data_local
                    wave_cnt[j] += len
                else
                    wave_data[j][:, wave_cnt[j] + (1:len - 1)] = wave_data_local[:, 2:end]
                    wave_cnt[j] += len - 1
                end

                if any(diff(wave_data[j][1, 1:wave_cnt[j]]) .<= 0.0) && wave_cnt[j] != length(unique(wave_data[j][1, 1:wave_cnt[j]]))
                    println("Warning: not all elements of the generated time vector are unique!")
                end
            end
        end
    end

    # trim the output data
    for j in 1:shape_channels
        if wave_cnt[j] < size(wave_data[j], 2)
            wave_data[j][:, (wave_cnt[j] + 1):end] = []
        end
    end

    return wave_data, tfp_excitation, tfp_refocusing, t_adc, fp_adc
end


function getBinaryCodes(obj)
    # Return binary codes for section headers in a binary sequence file.
    #
    # See also writeBinary

    codes = Dict{Symbol, Any}()
    codes[:fileHeader] = [1, "pulseq", 2]
    codes[:version_major] = Int64(obj.version_major)
    codes[:version_minor] = Int64(obj.version_minor)
    codes[:version_revision] = Int64(obj.version_revision)
    prefix = Int64(hex2dec("FFFFFFFFFFFFFFFF"))
    codes[:section] = Dict{Symbol, Int64}()
    codes[:section][:definitions] = prefix | Int64(1)
    codes[:section][:blocks] = prefix | Int64(2)
    codes[:section][:rf] = prefix | Int64(3)
    codes[:section][:gradients] = prefix | Int64(4)
    codes[:section][:trapezoids] = prefix | Int64(5)
    codes[:section][:adc] = prefix | Int64(6)
    codes[:section][:delays] = prefix | Int64(7)
    codes[:section][:shapes] = prefix | Int64(8)

    return codes
end

function getExtensionTypeID(obj, str)
    # get numeric ID for the given string extension ID
    # will automatically create a new ID if unknown

    num = findfirst(x -> x == str, obj.extensionStringIDs)

    if isnothing(num)
        if isempty(obj.extensionNumericIDs)
            id = 1
        else
            id = 1 + maximum(obj.extensionNumericIDs)
        end

        push!(obj.extensionNumericIDs, id)
        push!(obj.extensionStringIDs, str)
        @assert length(obj.extensionNumericIDs) == length(obj.extensionStringIDs)
    else
        id = obj.extensionNumericIDs[num]
    end

    return id
end

function getExtensionTypeString(obj, id)
    # get numeric ID for the given string extension ID
    # may fail

    num = findfirst(x -> x == id, obj.extensionNumericIDs)

    if isnothing(num)
        error("Extension for the given ID $id is unknown")
    end

    str = obj.extensionStringIDs[num]
    return str
end

function setExtensionStringAndID(obj, str, id)
    # set numeric ID for the given string extension ID
    # may fail if not unique

    if any(x -> x == str, obj.extensionStringIDs) || any(x -> x == id, obj.extensionNumericIDs)
        error("Numeric or String ID is not unique")
    end

    push!(obj.extensionNumericIDs, id)
    push!(obj.extensionStringIDs, str)

    @assert length(obj.extensionNumericIDs) == length(obj.extensionStringIDs)
end

#######################################################################
#######################################################################
#######################################################################

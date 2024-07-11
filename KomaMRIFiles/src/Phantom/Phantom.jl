"""
    phantom = read_phantom(filename)

Reads a (.phantom) file and creates a Phantom structure from it
"""
function read_phantom(filename::String)
    fid = HDF5.h5open(filename, "r")
    phantom_fields = []
    version = read_attribute(fid, "Version")
    dims = read_attribute(fid, "Dims")
    Ns = read_attribute(fid, "Ns")
    # Name 
    name = read_attribute(fid, "Name")
    push!(phantom_fields, (:name, name))
    # Position and contrast
    for key in ["position", "contrast"]
        group = fid[key]
        for label in HDF5.keys(group)
            param = group[label]
            values = read_param(param)
            if values != "Default"
                push!(phantom_fields, (Symbol(label), values))
            end
        end
    end
    # Motion
    motion_group = fid["motion"]
    import_motion!(phantom_fields, motion_group)

    obj = Phantom(; phantom_fields...)
    close(fid)
    return obj
end

function read_param(param::HDF5.Group)
    if "type" in HDF5.keys(attrs(param))
        type = attrs(param)["type"]
        if type == "Explicit"
            values = read(param["values"])
        elseif type == "Indexed"
            index = read(param["values"])
            @assert Ns == length(index) "Error: $(label) vector dimensions mismatch"
            table = read(param["table"])
            N = read_attribute(param, "N")
            @assert N == length(table) "Error: $(label) table dimensions mismatch"
            values = table[index]
        elseif type == "Default"
            values = "Default"
        end
    else
        values = read(param["values"])
    end
    return values
end

function import_motion!(phantom_fields::Array, motion_group::HDF5.Group)
    T = eltype(phantom_fields[2][2])
    motion_type = read_attribute(motion_group, "type")
    if motion_type == "MotionVector"
        simple_motion_types    = last.(split.(string.(reduce(vcat,(subtypes(subtypes(Motion)[2])))), "."))
        arbitrary_motion_types = last.(split.(string.(reduce(vcat,(subtypes(subtypes(Motion)[1])))), "."))
        motion_array = Motion{T}[]
        for key in keys(motion_group)
            type_group = motion_group[key]
            type_str = split(key, "_")[2]
            @assert type_str in vcat(simple_motion_types, arbitrary_motion_types) "Motion Type: $(type_str) has not been implemented in KomaMRIBase $(KomaMRIBase.__VERSION__)"
            args = []
            for smtype in subtypes(SimpleMotion)
                if type_str == last(split(string(smtype), "."))
                    times = import_time_range(type_group["times"])
                    type_fields = filter(x -> x != :times, fieldnames(smtype))
                    for key in type_fields
                        push!(args, read_attribute(type_group, string(key)))
                    end
                    push!(motion_array, smtype(times, args...))
                end
            end
            for amtype in subtypes(ArbitraryMotion)
                if type_str == last(split(string(amtype), "."))
                    times = import_time_range(type_group["times"])
                    type_fields = filter(x -> x != :times, fieldnames(amtype))
                    for key in type_fields
                        push!(args, read(type_group[string(key)]))
                    end
                    push!(motion_array, amtype(times, args...))
                end
            end
        end
        return push!(phantom_fields, (:motion, MotionVector(motion_array)))
    elseif motion_type == "NoMotion"
        return push!(phantom_fields, (:motion, NoMotion{T}()))
    end
end

function import_time_range(times_group::HDF5.Group)
    time_scale_type = read_attribute(times_group, "type")
    for tstype in subtypes(TimeScale)
        if time_scale_type == last(split(string(tstype), "."))
            args = []
            for key in filter(x -> x != :type, fieldnames(tstype))
                push!(args, read_attribute(times_group, string(key)))
            end
            return tstype(args...)
        end
    end
end

"""
    phantom = write_phantom(ph,filename)

Writes a (.phantom) file from a Phantom struct.
"""
function write_phantom(
    # By the moment, only "Explicit" type 
    # is considered when writing .phantom files
    obj::Phantom,
    filename::String;
    store_coords=[:x, :y, :z],
    store_contrasts=[:ρ, :T1, :T2, :T2s, :Δw],
    store_motion=true,
)
    # Create HDF5 phantom file
    fid = h5open(filename, "w")
    # Root attributes
    HDF5.attributes(fid)["Version"] = string(KomaMRIFiles.__VERSION__)
    HDF5.attributes(fid)["Name"] = obj.name
    HDF5.attributes(fid)["Ns"] = length(obj.x)
    dims = KomaMRIBase.get_dims(obj)
    HDF5.attributes(fid)["Dims"] = sum(dims)
    # Positions
    pos = create_group(fid, "position")
    for x in store_coords
        create_group(pos, String(x))["values"] = getfield(obj, x)
    end
    # Contrast (Rho, T1, T2, T2s Deltaw)
    contrast = create_group(fid, "contrast")
    for x in store_contrasts
        param = create_group(contrast, String(x))
        HDF5.attributes(param)["type"] = "Explicit" #TODO: consider "Indexed" type
        param["values"] = getfield(obj, x)
    end
    # Motion
    if store_motion
        motion_group = create_group(fid, "motion")
        export_motion!(motion_group, obj.motion)
    end
    return close(fid)
end

function export_motion!(motion_group::HDF5.Group, mv::MotionVector{T}) where {T<:Real}
    HDF5.attributes(motion_group)["type"] = "MotionVector"
    for (counter, m) in enumerate(mv.motions)
        type_name = typeof(m).name.name
        type_group = create_group(motion_group, "$(counter)_$type_name")
        export_time_range!(type_group, m.times)
        type_fields = filter(x -> x != :times, fieldnames(typeof(m)))
        for field in type_fields
            field_value = getfield(m, field)
            if typeof(field_value) <: Number
                HDF5.attributes(type_group)[string(field)] = field_value
            elseif typeof(field_value) <: AbstractArray
                type_group[string(field)] = field_value
            end
        end
    end
end

function export_motion!(motion_group::HDF5.Group, motion::NoMotion{T}) where {T<:Real}
    HDF5.attributes(motion_group)["type"] = "NoMotion"
end

function export_time_range!(type_group::HDF5.Group, times::TimeScale)
    times_name = typeof(times).name.name
    times_group = create_group(type_group, "times")
    HDF5.attributes(times_group)["type"] = string(times_name)
    for field in fieldnames(typeof(times))
        field_value = getfield(times, field)
        HDF5.attributes(times_group)[string(field)] = field_value
    end
end
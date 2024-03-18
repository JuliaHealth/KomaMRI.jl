"""
phantom = read_phantom(filename)

Reads a (.phantom) file and creates a Phantom structure from it
"""
function read_phantom(filename::String)
    # ----------------------------------------------
    function read_param(param::HDF5.Group)
        if "type" in HDF5.keys(attrs(param))
            type = attrs(param)["type"]

            if type == "Explicit"
                values = read(param["values"])
            elseif type == "Indexed"
                index = read(param["values"])
                if Ns == length(index)
                    table = read(param["table"])
                    N = read_attribute(param, "N")
                    if N == length(table)
                        values = table[index]
                    else
                        print("Error: $(label) table dimensions mismatch")
                    end
                else
                    print("Error: $(label) vector dimensions mismatch")
                end
            elseif type == "Default"
                values = "Default"
            end
        else
            values = read(param["values"])
        end

        return values
    end
    # ----------------------------------------------

    fid = HDF5.h5open(filename, "r")

    name = read_attribute(fid, "Name")
    version = read_attribute(fid, "Version")
    dims = read_attribute(fid, "Dims")
    Ns = read_attribute(fid, "Ns")
    precision = read_attribute(fid, "Precision")
    if precision == "f32"
        tp = Float32
    elseif precision == "f64"
        tp = Float64
    end
    obj = Phantom{tp}(; name=name, x=zeros(Ns))

    # Position and contrast
    for key in ["position", "contrast"]
        group = fid[key]
        for label in HDF5.keys(group)
            param = group[label]
            values = read_param(param)
            if values != "Default"
                setfield!(obj, Symbol(label), values)
            end
        end
    end

    # Diffusion (TODO)

    # Motion
    motion_group = fid["motion"]
    obj.motion = import_motion(motion_group, tp)

    close(fid)
    return obj
end

function import_motion(motion_group::HDF5.Group, precision::Type)
    model = read_attribute(motion_group, "model")
    if model == "NoMotion"
        return NoMotion()
    elseif model == "SimpleMotion"
        motion = SimpleMotion(SimpleMotionType{precision}[])
        types_group = motion_group["types"]
        for key in keys(types_group)
            type_group = types_group[key]
            type_str = match(r"^\d+_(\w+)", key).captures[1]
            type = type_str == "Translation" ? Translation :
                   type_str == "Rotation" ? Rotation : 
                   nothing
            args = []
            for key in fieldnames(type)
                push!(args, read_attribute(type_group, string(key)))
            end
            
            push!(motion.types, type(args...))
        end
        return motion

    #TODO: import_motion for ArbitraryMotion
    elseif model == "ArbitraryMotion"
    #     segments = motion["segments"]
    #     N = read_attribute(segments, "N")
    #     K = read_attribute(segments, "K")
    #     dur = read(segments["dur"])

    #     Δx = zeros(Ns, K - 1)
    #     Δy = zeros(Ns, K - 1)
    #     Δz = zeros(Ns, K - 1)

    #     resetmag = zeros(Ns, K)

    #     for key in HDF5.keys(motion)
    #         if key != "segments"
    #             values = read_param(motion[key])
    #             if key == "motionx"
    #                 Δx = values
    #             elseif key == "motiony"
    #                 Δy = values
    #             elseif key == "motionz"
    #                 Δz = values
    #             end
    #         end
    #     end
    end
end



"""
write_phantom(ph,filename)

Writes a (.phantom) file from a Phantom struct.
"""
# By the moment, only "Explicit" type 
# is considered when writing .phantom files
function write_phantom(obj::Phantom, filename::String)
    # Create HDF5 phantom file
    fid = h5open(filename, "w")

    # Root attributes
    HDF5.attributes(fid)["Version"] = "0.1"
    HDF5.attributes(fid)["Name"] = obj.name
    HDF5.attributes(fid)["Ns"] = length(obj.x)
    dims = get_dims(obj)
    HDF5.attributes(fid)["Dims"] = sum(dims)
    # Precision
    if eltype(obj.x) <: Float32
        HDF5.attributes(fid)["Precision"] = "f32"
    elseif eltype(obj.x) <: Float64
        HDF5.attributes(fid)["Precision"] = "f64"
    end

    fields = fieldnames(Phantom)[2:end]

    # Spin initial positions
    pos = create_group(fid, "position")
    for i in 1:3
        if Bool(dims[i])
            create_group(pos, String(fields[i]))["values"] = getfield(obj, fields[i])
        end
    end

    # Contrast (Rho, T1, T2, T2s Deltaw)
    contrast = create_group(fid, "contrast")
    for i in 4:8
        param = create_group(contrast, String(fields[i]))
        HDF5.attributes(param)["type"] = "Explicit" #TODO: consider "Indexed" type
        param["values"] = getfield(obj, fields[i])
    end

    # Diffusion (TODO)

    # Motion
    motion = create_group(fid, "motion")
    export_motion(motion, obj.motion)

    return close(fid)
end

function export_motion(motion_group::HDF5.Group, motion::NoMotion)
    HDF5.attributes(motion_group)["model"] = "NoMotion"
end

function export_motion(motion_group::HDF5.Group, motion::SimpleMotion)
    HDF5.attributes(motion_group)["model"] = "SimpleMotion"
    types_group =  create_group(motion_group, "types")
    counter = 1
    for type in motion.types
        type_group = create_group(types_group, string(counter)*"_"*match(r"^\w+", string(typeof(type))).match)
        fields = fieldnames(typeof(type))
        for i in 1:length(fields)
            HDF5.attributes(type_group)[string(fields[i])] = getfield(type, fields[i])
        end
        counter += 1
    end
end

#TODO: export_motion for ArbitraryMotion
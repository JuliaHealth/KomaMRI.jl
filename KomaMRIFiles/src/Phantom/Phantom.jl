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
    model = read_attribute(motion_group, "model")
    import_motion!(phantom_fields, Ns, Symbol(model), motion_group)

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

function import_motion!(phantom_fields::Array, Ns::Int, model::Symbol, motion_group::HDF5.Group)
    return import_motion!(phantom_fields, Ns, Val(model), motion_group)
end
function import_motion!(
    phantom_fields::Array, Ns::Int, model::Val{:NoMotion}, motion_group::HDF5.Group
)
    return nothing
end
function import_motion!(
    phantom_fields::Array, Ns::Int, model::Val{:SimpleMotion}, motion_group::HDF5.Group
)
    types_group = motion_group["types"]
    types = SimpleMotionType[]
    for key in keys(types_group)
        type_group = types_group[key]
        type_str = split(key, "_")[2]
        @assert type_str in last.(split.(string.(subtypes(SimpleMotionType)), ".")) "Simple Motion Type: $(type_str) has not been implemented in KomaMRIBase $(KomaMRIBase.__VERSION__)"
        for SMT in subtypes(SimpleMotionType)
            args = []
            if type_str == last(split(string(SMT), "."))
                for key in fieldnames(SMT)
                    push!(args, read_attribute(type_group, string(key)))
                end
                push!(types, SMT(args...))
            end
        end
    end
    return push!(phantom_fields, (:motion, SimpleMotion((types...))))
end
function import_motion!(
    phantom_fields::Array, Ns::Int, model::Val{:ArbitraryMotion}, motion_group::HDF5.Group
)
    t_start = read(motion_group["t_start"])
    t_end = read(motion_group["t_end"])
    dx = read(motion_group["dx"])  
    dy = read(motion_group["dy"])  
    dz = read(motion_group["dz"])  
    return push!(phantom_fields, (:motion, ArbitraryMotion(t_start, t_end, dx, dy, dz)))
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

function export_motion!(motion_group::HDF5.Group, motion::NoMotion)
    return HDF5.attributes(motion_group)["model"] = "NoMotion"
end
function export_motion!(motion_group::HDF5.Group, motion::SimpleMotion)
    HDF5.attributes(motion_group)["model"] = "SimpleMotion"
    types_group = create_group(motion_group, "types")
    counter = 1
    for (counter, sm_type) in enumerate(motion.types)
        simple_motion_type = typeof(sm_type).name.name
        type_group = create_group(types_group, "$(counter)_$simple_motion_type")
        phantom_fields = fieldnames(typeof(sm_type))
        for field in phantom_fields
            HDF5.attributes(type_group)[string(field)] = getfield(sm_type, field)
        end
    end
end
function export_motion!(motion_group::HDF5.Group, motion::ArbitraryMotion)
    HDF5.attributes(motion_group)["model"] = "ArbitraryMotion"
    motion_group["t_start"] = motion.t_start
    motion_group["t_end"] = motion.t_end
    motion_group["dx"] = motion.dx
    motion_group["dy"] = motion.dy
    motion_group["dz"] = motion.dz
end
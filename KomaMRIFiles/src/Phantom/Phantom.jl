"""
    phantom = read_phantom(filename)

Reads a (.phantom) file and creates a Phantom structure from it
"""
function read_phantom(filename::String)
    fid = HDF5.h5open(filename, "r")
    dims = read_attribute(fid, "Dims")
    Ns = read_attribute(fid, "Ns")
    # Version
    file_version = VersionNumber(read_attribute(fid, "Version"))
    program_version = pkgversion(KomaMRIFiles)
    if file_version.major != program_version.major
        @warn "KomaMRIFiles: Version mismatch detected: $file_version (used to write .phantom) vs $program_version (installed)
         This may lead to compatibility issues. "
    end
    phantom_fields = []
    # Name 
    name = read_attribute(fid, "Name")
    push!(phantom_fields, (:name, name))
    # Position and contrast
    for key in ["position", "contrast"]
        group = fid[key]
        for label in HDF5.keys(group)
            values = read(group[label])
            push!(phantom_fields, (Symbol(label), values))
        end
    end
    # Motion
    if "motion" in keys(fid)
        motion_group = fid["motion"]
        import_motion!(phantom_fields, motion_group)
    end
    obj = Phantom(; phantom_fields...)
    close(fid)
    return obj
end

function import_motion!(phantom_fields::Array, motion_group::HDF5.Group)
    T = eltype(phantom_fields[2][2])
    motion_array = Motion{T}[]
    for key in keys(motion_group)
        motion = motion_group[key]
        motion_fields = []
        for name in keys(motion) # action, time, spins
            import_motion_field!(motion_fields, motion, name)
        end
        push!(motion_array, Motion(; motion_fields...))
    end
    push!(phantom_fields, (:motion, MotionList(motion_array)))
end

function import_motion_field!(motion_fields::Array, motion::HDF5.Group, name::String)
    field_group = motion[name]
    type = read_attribute(field_group, "type")

    get_subtypes(t::Type) = reduce(vcat,(subtypes(t)))
    get_subtype_strings(t::Type) = last.(split.(string.(get_subtypes(t::Type)), "."))
    
    subtype_strings = reduce(vcat, get_subtype_strings.([
        KomaMRIBase.SimpleAction,
        KomaMRIBase.ArbitraryAction,
        KomaMRIBase.AbstractTimeSpan,
        KomaMRIBase.AbstractSpinSpan
    ]))

    subtype_vector = reduce(vcat, get_subtypes.([
        KomaMRIBase.SimpleAction,
        KomaMRIBase.ArbitraryAction,
        KomaMRIBase.AbstractTimeSpan,
        KomaMRIBase.AbstractSpinSpan
    ]))

    motion_subfields = []
    for (i, subtype_string) in enumerate(subtype_strings)
        if type == subtype_string
            for subname in fieldnames(subtype_vector[i]) # dx, dy, dz, pitch, roll...
                key = string(subname)
                subfield_value = key in keys(field_group) ? read(field_group, key) : read_attribute(field_group, key)
                import_motion_subfield!(motion_subfields, subfield_value, key)
            end
            push!(motion_fields, (Symbol(name), subtype_vector[i](motion_subfields...)))
        end
    end
end

function import_motion_subfield!(motion_subfields::Array, subfield_value::Union{Real, Array}, key::String)
    push!(motion_subfields, subfield_value)
    return nothing
end
function import_motion_subfield!(motion_subfields::Array, subfield_value::String, key::String)
    endpoints = parse.(Int, split(subfield_value, ":"))
    range = length(endpoints) == 3 ? (endpoints[1]:endpoints[2]:endpoints[3]) : (endpoints[1]:endpoints[2])
    push!(motion_subfields, range)
    return nothing
end


"""
    phantom = write_phantom(ph,filename)

Writes a (.phantom) file from a Phantom struct.
"""
function write_phantom(
    obj::Phantom,
    filename::String;
    store_coords=[:x, :y, :z],
    store_contrasts=[:ρ, :T1, :T2, :T2s, :Δw],
    store_motion=true
)
    # Create HDF5 phantom file
    fid = h5open(filename, "w")
    # Root attributes
    HDF5.attributes(fid)["Version"] = string(pkgversion(KomaMRIFiles))
    HDF5.attributes(fid)["Name"] = obj.name
    HDF5.attributes(fid)["Ns"] = length(obj.x)
    dims = KomaMRIBase.get_dims(obj)
    HDF5.attributes(fid)["Dims"] = sum(dims)
    # Positions
    pos = create_group(fid, "position")
    for x in store_coords
        pos[String(x)] = getfield(obj, x)
    end
    # Contrast (Rho, T1, T2, T2s Deltaw)
    contrast = create_group(fid, "contrast")
    for x in store_contrasts
        contrast[String(x)] = getfield(obj, x)
    end
    # Motion
    if (typeof(obj.motion) <: MotionList) & store_motion
        motion_group = create_group(fid, "motion")
        export_motion!(motion_group, obj.motion)
    end
    return close(fid)
end

function export_motion!(motion_group::HDF5.Group, motion_list::MotionList)
    KomaMRIBase.sort_motions!(motion_list)
    for (counter, m) in enumerate(motion_list.motions)
        motion = create_group(motion_group, "motion_$(counter)")
        for key in fieldnames(Motion) # action, time, spins
            field_group = create_group(motion, string(key)) 
            field_value = getfield(m, key)
            export_motion_field!(field_group, field_value)
        end
    end
end

function export_motion_field!(field_group::HDF5.Group, field_value)
    HDF5.attributes(field_group)["type"] = string(typeof(field_value).name.name)
    for subname in fieldnames(typeof(field_value)) # dx, dy, dz, pitch, roll...
        subfield = getfield(field_value, subname)
        export_motion_subfield!(field_group, subfield, string(subname))
    end
end

function export_motion_subfield!(field_group::HDF5.Group, subfield::Real, subname::String)
    HDF5.attributes(field_group)[subname] = subfield
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::AbstractRange, subname::String)
    HDF5.attributes(field_group)[subname] = step(subfield) == 1 ? "$(first(subfield)):$(last(subfield))" : "$(first(subfield)):$(step(subfield)):$(last(subfield))"
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::Array, subname::String)
    field_group[subname] = subfield
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::BitMatrix, subname::String)
    field_group[subname] = Int.(subfield)
end
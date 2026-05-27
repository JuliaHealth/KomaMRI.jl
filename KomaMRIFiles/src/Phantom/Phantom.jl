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
    pos = fid["position"]
    for label in HDF5.keys(pos)
        push!(phantom_fields, (Symbol(label), read(pos[label])))
    end
    contrast = fid["contrast"]
    T = eltype(phantom_fields[end][2])
    for label in HDF5.keys(contrast)
        import_contrast_property!(phantom_fields, contrast, label, T)
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
            import_motion_field!(motion_fields, motion, name, T)
        end
        push!(motion_array, Motion(; motion_fields...))
    end
    push!(phantom_fields, (:motion, MotionList(motion_array...)))
end

function import_motion_field!(motion_fields::Array, motion::HDF5.Group, name::String, T::Type{<:Real})
    field_group = motion[name]
    type = read_attribute(field_group, "type")

    get_subtypes(t::Type) = reduce(vcat,(subtypes(t)))
    get_subtype_strings(t::Type) = last.(split.(string.(get_subtypes(t::Type)), "."))
    
    subtype_strings = [reduce(vcat, get_subtype_strings.([
        SimpleAction,
        ArbitraryAction,
        AbstractSpinSpan
    ])); "TimeCurve"] 

    subtype_vector = [reduce(vcat, get_subtypes.([
        SimpleAction,
        ArbitraryAction,
        AbstractSpinSpan
    ])); TimeCurve]

    motion_subfields = []
    for (i, subtype_string) in enumerate(subtype_strings)
        if type == subtype_string
            for subname in fieldnames(subtype_vector[i]) # dx, dy, dz, pitch, roll...
                key = string(subname)
                if !(key in ["t_start", "t_end"])
                    subfield_value = key in keys(field_group) ? read(field_group, key) : read_attribute(field_group, key)
                    import_motion_subfield!(motion_subfields, subfield_value, key, T)
                end
            end
            push!(motion_fields, (Symbol(name), subtype_vector[i](motion_subfields...)))
        end
    end
end

function import_motion_subfield!(motion_subfields::Array, subfield_value::Union{Real, Array}, key::String, T::Type{<:Real})
    push!(motion_subfields, subfield_value)
    return nothing
end
function read_timecurve_group(group::HDF5.Group)
    t = read(group, "t")
    t_unit = read(group, "t_unit")
    periodic = haskey(HDF5.attributes(group), "periodic") ? read_attribute(group, "periodic") : false
    periods = haskey(HDF5.attributes(group), "periods") ? read_attribute(group, "periods") : 1.0
    return TimeCurve(; t, t_unit, periodic, periods)
end

function write_timecurve_group!(group::HDF5.Group, tc::TimeCurve)
    HDF5.attributes(group)["type"] = "TimeCurve"
    group["t"] = tc.t
    group["t_unit"] = tc.t_unit
    HDF5.attributes(group)["periodic"] = tc.periodic
    HDF5.attributes(group)["periods"] = tc.periods
end

function export_contrast_property!(contrast::HDF5.Group, label::Symbol, prop)
    if prop isa TimeDependentProperty
        grp = create_group(contrast, string(label))
        HDF5.attributes(grp)["type"] = "TimeDependentProperty"
        grp["value"] = prop.value
        time_grp = create_group(grp, "time")
        write_timecurve_group!(time_grp, prop.time)
    else
        contrast[string(label)] = prop
    end
    return nothing
end

function import_contrast_property!(phantom_fields::Array, contrast::HDF5.Group, label::String, ::Type{T}) where {T<:Real}
    node = contrast[label]
    if node isa HDF5.Dataset
        push!(phantom_fields, (Symbol(label), read(node)))
    else
        @assert read_attribute(node, "type") == "TimeDependentProperty"
        value = read(node["value"])
        time = read_timecurve_group(node["time"])
        push!(phantom_fields, (Symbol(label), TimeDependentProperty(value, time)))
    end
    return nothing
end

function import_motion_subfield!(motion_subfields::Array, subfield_value::String, key::String, T::Type{<:Real})
    if subfield_value in ["true", "false"]
        return push!(motion_subfields, subfield_value == "true" ? true : false)
    elseif subfield_value == "CenterOfMass"
        return push!(motion_subfields, CenterOfMass())
    elseif startswith(subfield_value, "(") && endswith(subfield_value, ")")
        elements = split(subfield_value[2:end-1], ",")
        return push!(motion_subfields, Tuple(parse.(T, strip.(elements))))
    end
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
    dims = get_dims(obj)
    HDF5.attributes(fid)["Dims"] = sum(dims)
    # Positions
    pos = create_group(fid, "position")
    for x in store_coords
        pos[String(x)] = getfield(obj, x)
    end
    # Contrast (Rho, T1, T2, T2s Deltaw)
    contrast = create_group(fid, "contrast")
    for x in store_contrasts
        export_contrast_property!(contrast, x, getfield(obj, x))
    end
    # Motion
    if !(obj.motion isa NoMotion) & store_motion
        motion_group = create_group(fid, "motion")
        export_motion!(motion_group, obj.motion)
    end
    return close(fid)
end

export_motion!(motion_group::HDF5.Group, motion::Motion) = export_motion!(motion_group, MotionList([motion]))
function export_motion!(motion_group::HDF5.Group, motion_list::MotionList)
    sort_motions!(motion_list)
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
        !(string(subname) in ["t_start", "t_end"]) ? export_motion_subfield!(field_group, subfield, string(subname)) : nothing
    end
end

function export_motion_subfield!(field_group::HDF5.Group, subfield::Real, subname::String)
    HDF5.attributes(field_group)[subname] = subfield
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::AbstractRange, subname::String)
    HDF5.attributes(field_group)[subname] = step(subfield) == 1 ? "$(first(subfield)):$(last(subfield))" : "$(first(subfield)):$(step(subfield)):$(last(subfield))"
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::Union{Bool,Tuple}, subname::String)
    HDF5.attributes(field_group)[subname] = string(subfield)
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::Array, subname::String)
    field_group[subname] = subfield
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::BitMatrix, subname::String)
    field_group[subname] = Int.(subfield)
end
function export_motion_subfield!(field_group::HDF5.Group, subfield::CenterOfMass, subname::String)
    field_group[subname] = "CenterOfMass"
end

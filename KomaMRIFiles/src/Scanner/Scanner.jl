"""
    write_scanner(sys, filename)

Write a scanner to an HDF5 `.sys` file. The `limits`, `gradient`, `receiver`, and
`transmitter` groups mirror the fields of `Scanner`. Each group has a `type`
attribute, scalar parameters as attributes, and arrays as datasets. Coordinates
and hardware limits use SI units; complex receive maps retain their values and
`(x, y, z, coil)` dimensions.

Supports `HardwareLimits`, `LinearXYZ`, `UniformCoilSens`, `BirdcageCoilSens`,
`ArbitraryCoilSens`, and `UniformTransmit`. The root `Version` attribute records
the KomaMRIFiles version, as in `.phantom` files.
"""
function write_scanner(sys, filename)
    h5open(filename, "w") do file
        attributes(file)["Version"] = string(pkgversion(KomaMRIFiles))
        for name in fieldnames(Scanner)
            group = create_group(file, string(name))
            write_scanner_component!(group, getproperty(sys, name))
        end
    end
    return nothing
end

function write_scanner_component!(group, component)
    attributes(group)["type"] = string(nameof(typeof(component)))
    for name in fieldnames(typeof(component))
        write_scanner_field!(group, string(name), getproperty(component, name))
    end
    return nothing
end

function write_scanner_field!(group, name, value::Number)
    attributes(group)[name] = value
    return nothing
end

function write_scanner_field!(group, name, value::AbstractArray)
    group[name] = Array(value)
    return nothing
end

"""
    sys = read_scanner(filename)

Read an HDF5 `.sys` file written by [`write_scanner`](@ref), reconstructing the
scanner's hardware limits and gradient, receive, and transmit models. Model names
are resolved from the loaded scanner subtypes, as for phantom motion; file
contents are never evaluated as Julia code.
"""
function read_scanner(filename)
    return h5open(filename, "r") do file
        version = VersionNumber(read_attribute(file, "Version"))
        version.major == pkgversion(KomaMRIFiles).major ||
            @warn "Scanner file was written with a different KomaMRIFiles major version" file_version=version installed_version=pkgversion(KomaMRIFiles)
        components = map(fieldnames(Scanner)) do name
            name => read_scanner_component(file[string(name)])
        end
        Scanner(; components...)
    end
end

function read_scanner_component(group)
    name = read_attribute(group, "type")
    types = (HardwareLimits, subtypes(AbstractGradientSystem)...,
        subtypes(AbstractRFReceiveSystem)..., subtypes(AbstractRFTransmitSystem)...)
    index = findfirst(type -> string(nameof(type)) == name, types)
    isnothing(index) && throw(ArgumentError("Unsupported scanner component: $name"))
    type = types[index]
    values = map(fieldnames(type)) do field
        key = string(field)
        haskey(group, key) ? read(group[key]) : read_attribute(group, key)
    end
    return type(values...)
end

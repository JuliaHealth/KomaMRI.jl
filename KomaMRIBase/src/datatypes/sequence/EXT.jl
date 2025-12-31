include("extensions/AdcLabels.jl")
export AdcLabels

abstract type Extension end
# Supported extensions. To add a new extension, create a new file in the extensions folder and add it to the list below.
include("extensions/LabelInc.jl")
include("extensions/LabelSet.jl")
include("extensions/Trigger.jl")

get_EXT_type_from_symbol(::Val{:LABELINC}) = LabelInc
get_EXT_type_from_symbol(::Val{:LABELSET}) = LabelSet
get_EXT_type_from_symbol(::Val{:TRIGGERS}) = Trigger
get_EXT_type_from_symbol(::Val) = nothing

"""
    format_string = get_scanf_format(T)

Generates a scanf format string from the field types of struct type `T`.
The format string is generated in the order fields are defined in the struct.

"""
function get_scanf_format(::Type{T}) where T
    type_to_format = Dict(
        Int => "%i",
        Int64 => "%i",
        Int32 => "%i",
        Int16 => "%i",
        Int8 => "%i",
        Float64 => "%f",
        Float32 => "%f",
        String => "%s",
        Char => "%c"
    )
    format_parts = [get(type_to_format, ft, "%s") for ft in fieldtypes(T)]
    return join(format_parts, " ")
end

export Extension, LabelInc, LabelSet, Trigger, SoftDelay
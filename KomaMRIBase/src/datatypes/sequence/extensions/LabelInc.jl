mutable struct LabelInc <: Extension
    labelvalue::Int
    labelstring::String
end

get_scale(::Type{LabelInc}) = [1 1]
get_scanf_format(::Type{LabelInc}) = "%i %s"
get_EXT_type_from_symbol(::Val{:LABELINC}) = LabelInc
extension_type_header(::Type{LabelInc}) = "# Extension specification for setting labels:\n# id set labelstring\n"
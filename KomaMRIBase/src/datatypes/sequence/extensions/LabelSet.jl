mutable struct LabelSet <: Extension
    labelvalue::Int
    labelstring::String
end

get_scale(::Type{LabelSet}) = [1 1]
get_scanf_format(::Type{LabelSet}) = "%i %s"
get_EXT_type_from_symbol(::Val{:LABELSET}) = LabelSet
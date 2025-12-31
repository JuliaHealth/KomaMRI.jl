mutable struct LabelSet <: Extension
    labelvalue::Int
    labelstring::String
end

get_scale(::Type{LabelSet}) = [1 ""]
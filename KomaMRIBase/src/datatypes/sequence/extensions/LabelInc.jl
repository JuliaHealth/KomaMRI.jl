mutable struct LabelInc <: Extension
    labelvalue::Int
    labelstring::String
end

get_scale(::Type{LabelInc}) = [1 ""]
struct LabelInc <: Extension
    labelvalue::Int
    labelstring::String
end

get_scale(::Type{LabelInc}) = [1 1]
get_pulseq_format(::Type{LabelInc}) = "%i %s"
get_EXT_type_from_symbol(::Val{:LABELINC}) = LabelInc
get_symbol_from_EXT_type(::Type{LabelInc}) = "LABELINC"
extension_type_header(::Type{LabelInc}) = "# Extension specification for setting labels:\n# id set labelstring\n"

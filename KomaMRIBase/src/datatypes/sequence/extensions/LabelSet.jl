struct LabelSet <: Extension
    labelvalue::Int
    labelstring::String
end

get_scale(::Type{LabelSet}) = [1 1]
get_pulseq_format(::Type{LabelSet}) = "%i %s"
get_EXT_type_from_symbol(::Val{:LABELSET}) = LabelSet
get_symbol_from_EXT_type(::Type{LabelSet}) = "LABELSET"
extension_type_header(::Type{LabelSet}) = "# Extension specification for setting labels:\n# id set labelstring\n"

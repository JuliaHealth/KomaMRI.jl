include("extensions/AdcLabels.jl")
export AdcLabels

abstract type Extension end
# Supported extensions. To add a new extension, create a new file in the extensions folder and add it to the list below.
include("extensions/LabelInc.jl")
include("extensions/LabelSet.jl")
include("extensions/Trigger.jl")

get_EXT_type_from_symbol(::Val) = nothing

export Extension, LabelInc, LabelSet, Trigger, SoftDelay
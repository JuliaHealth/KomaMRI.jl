include("extensions/AdcLabels.jl")
export AdcLabels

abstract type Extension end

# Supported extensions. To add a new extension, create a new file in the extensions folder and add it below.
include("extensions/LabelInc.jl")
include("extensions/LabelSet.jl")
include("extensions/Trigger.jl")
include("extensions/QuaternionRot.jl")

get_EXT_type_from_symbol(::Val) = nothing

export Extension, LabelInc, LabelSet, Trigger, SoftDelay, QuaternionRot, apply_rotations!

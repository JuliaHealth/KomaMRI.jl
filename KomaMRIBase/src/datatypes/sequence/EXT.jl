include("extensions/AdcLabels.jl")
export AdcLabels

abstract type Extension end
# Supported extensions. To add a new extension, create a new file in the extensions folder and add it to the list below.
include("extensions/LabelInc.jl")
include("extensions/LabelSet.jl")
include("extensions/Trigger.jl")

get_EXT_type_from_symbol(::Val) = nothing

Base.:(==)(x::Extension, y::Extension) = (typeof(x) == typeof(y)) && all([getfield(x, k) == getfield(y, k) for k in fieldnames(typeof(x))])
function Base.:(≈)(x::Extension, y::Extension)
    typeof(x) == typeof(y) || return false
    equal = true
    for k in fieldnames(typeof(x))
        fx, fy = getfield(x, k), getfield(y, k)
        equal &= fx isa Number ? fx ≈ fy : fx == fy
    end
    return equal
end

export Extension, LabelInc, LabelSet, Trigger, SoftDelay
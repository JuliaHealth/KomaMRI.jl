include("extensions/AdcLabels.jl")
export AdcLabels

abstract type Extension end
dur(::Extension) = 0.0

# Supported extensions. To add a new extension, create a new file in the extensions folder and add it to the list below.
include("extensions/LabelInc.jl")
include("extensions/LabelSet.jl")
include("extensions/QuaternionRot.jl")
include("extensions/Trigger.jl")

get_EXT_type_from_symbol(::Val)  = nothing
get_symbol_from_EXT_type(::Type) = nothing

Base.:(==)(::Extension, ::Extension) = false
Base.:(==)(x::T, y::T) where {T<:Extension} = fields_equal(x, y)

Base.isapprox(::Extension, ::Extension; kwargs...) = false
field_isapprox(x::Extension, y::Extension; kwargs...) = isapprox(x, y; kwargs...)
Base.isapprox(x::T, y::T; kwargs...) where {T<:Extension} = fields_isapprox(x, y; kwargs...)

export Extension, LabelInc, LabelSet, QuaternionRot, Trigger

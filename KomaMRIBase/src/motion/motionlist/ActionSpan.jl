abstract type AbstractActionSpan{T<:Real} end

is_composable(m::AbstractActionSpan) = false

# Simple actions
include("actions/SimpleAction.jl")
# Arbitrary actions
include("actions/ArbitraryAction.jl")

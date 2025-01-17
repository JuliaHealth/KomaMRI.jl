abstract type AbstractAction{T<:Real} end

is_composable(m::AbstractAction) = true

# Simple actions
include("actions/SimpleAction.jl")
# Arbitrary actions
include("actions/ArbitraryAction.jl")

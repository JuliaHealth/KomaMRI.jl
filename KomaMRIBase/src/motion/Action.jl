abstract type AbstractAction{T<:Real} end

Base.:(≈)(a, b) = (a isa Nothing || b isa Nothing) ? false : Base.≈(a, b)
Base.:(≈)(::Nothing, ::Nothing) = true
Base.:(≈)(a1::AbstractAction, a2::AbstractAction) = (typeof(a1) == typeof(a2)) & reduce(&, [getfield(a1, field)  ≈ getfield(a2, field) for field in fieldnames(typeof(a1))])
is_composable(m::AbstractAction) = true

# Simple actions
include("actions/SimpleAction.jl")
# Arbitrary actions
include("actions/ArbitraryAction.jl")

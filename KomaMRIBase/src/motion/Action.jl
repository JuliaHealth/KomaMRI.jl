abstract type AbstractAction{T<:Real} end

Base.:(==)(a1::AbstractAction, a2::AbstractAction) = (typeof(a1) == typeof(a2)) & reduce(&, [getfield(a1, field) == getfield(a2, field) for field in fieldnames(typeof(a1))])
Base.:(≈)(a1::AbstractAction,  a2::AbstractAction) = (typeof(a1) == typeof(a2)) & reduce(&, [getfield(a1, field)  ≈ getfield(a2, field) for field in fieldnames(typeof(a1))])
is_composable(::AbstractAction) = true
has_cycle_map(::AbstractAction) = false
add_reset_times!(t, ::AbstractAction, t_start, t_end, periods) = nothing
add_cycle_remap_times!(t, ::AbstractAction, t_start, t_end, periods) = nothing
cycle_remap_sources(::AbstractAction, spins, x) = nothing

# Simple actions
include("actions/SimpleAction.jl")
# Arbitrary actions
include("actions/ArbitraryAction.jl")

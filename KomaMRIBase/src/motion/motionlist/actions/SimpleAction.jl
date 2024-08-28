abstract type SimpleAction{T<:Real} <: AbstractActionSpan{T} end

Base.getindex(action::SimpleAction, p::Union{AbstractVector, Colon}) = action
Base.view(action::SimpleAction, p::Union{AbstractVector, Colon}) = action

include("simpleactions/Translate.jl")
include("simpleactions/Rotate.jl")
include("simpleactions/HeartBeat.jl")
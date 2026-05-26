abstract type SimpleAction <: AbstractAction end

Base.getindex(action::SimpleAction, p) = action
Base.view(action::SimpleAction, p)     = action

include("simpleactions/Translate.jl")
include("simpleactions/Rotate.jl")
include("simpleactions/HeartBeat.jl")
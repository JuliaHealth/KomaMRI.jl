abstract type ArbitraryAction{T<:Real} <: AbstractAction{T} end

function Base.getindex(action::ArbitraryAction, p)
    return typeof(action)([getfield(action, d)[p,:] for d in fieldnames(typeof(action))]...)
end
function Base.view(action::ArbitraryAction, p)
    return typeof(action)([@view(getfield(action, d)[p,:]) for d in fieldnames(typeof(action))]...)
end

function displacement_x!(ux, action::ArbitraryAction, x, y, z, t)
    resample_linear!(ux, action.dx, t)
    return nothing
end

function displacement_y!(uy, action::ArbitraryAction, x, y, z, t)
    resample_linear!(uy, action.dy, t)
    return nothing
end

function displacement_z!(uz, action::ArbitraryAction, x, y, z, t)
    resample_linear!(uz, action.dz, t)
    return nothing
end

include("arbitraryactions/Path.jl")
include("arbitraryactions/FlowPath.jl")

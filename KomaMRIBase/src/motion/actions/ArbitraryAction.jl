abstract type ArbitraryAction{T<:Real} <: AbstractAction{T} end

function Base.getindex(action::ArbitraryAction, p)
    return typeof(action)([getfield(action, d)[p,:] for d in fieldnames(typeof(action))]...)
end
function Base.view(action::ArbitraryAction, p)
    return typeof(action)([@view(getfield(action, d)[p,:]) for d in fieldnames(typeof(action))]...)
end

function displacement_x!(ux, action::ArbitraryAction, x, y, z, t)
    itp = interpolate(action.dx, Gridded(Linear()), Val(size(action.dx,1)), t)
    ux .= resample(itp, t)
    return nothing
end

function displacement_y!(uy, action::ArbitraryAction, x, y, z, t)
    itp = interpolate(action.dy, Gridded(Linear()), Val(size(action.dy,1)), t)
    uy .= resample(itp, t)
    return nothing
end

function displacement_z!(uz, action::ArbitraryAction, x, y, z, t)
    itp = interpolate(action.dz, Gridded(Linear()), Val(size(action.dz,1)), t)
    uz .= resample(itp, t)
    return nothing
end

include("arbitraryactions/Path.jl")
include("arbitraryactions/FlowPath.jl")
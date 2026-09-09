"""Receive model with one spatially uniform channel."""
struct UniformCoilSens <: AbstractRFReceiveSystem end

get_n_coils(::UniformCoilSens) = 1

function get_sens(receiver::UniformCoilSens, x, y, z)
    sens = similar(x, Complex{eltype(x)}, length(x), 1)
    return get_sens!(sens, receiver, x, y, z)
end

function get_sens!(sens, ::UniformCoilSens, x, y, z)
    fill!(sens, one(eltype(sens)))
    return sens
end

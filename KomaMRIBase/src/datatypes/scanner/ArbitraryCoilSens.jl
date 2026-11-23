"""
    ArbitraryCoilSens(x, y, z, coil_sens)

Receive model defined by complex sensitivity samples on a Cartesian grid.
`coil_sens` must have size `(length(x), length(y), length(z), ncoils)`.
Non-singleton coordinate vectors must be strictly increasing.
Sensitivities outside the sampled grid are zero.
"""
struct ArbitraryCoilSens{
    X<:AbstractVector,
    Y<:AbstractVector,
    Z<:AbstractVector,
    S<:AbstractArray{<:Complex,4},
} <: AbstractRFReceiveSystem
    x::X
    y::Y
    z::Z
    coil_sens::S
end

get_n_coils(receiver::ArbitraryCoilSens) = size(receiver.coil_sens, 4)

function coil_interpolators(receiver::ArbitraryCoilSens)
    grids = (receiver.x, receiver.y, receiver.z)
    interpolated_dims = findall(>(1), length.(grids))
    singleton_dims = findall(==(1), length.(grids))

    isempty(interpolated_dims) && return (; interpolated_dims, singleton_dims)
    knots = Tuple(grids[dim] for dim in interpolated_dims)
    interpolations = [
        linear_interpolation(
            knots,
            dropdims(
                @view(receiver.coil_sens[:, :, :, coil]); dims=Tuple(singleton_dims),
            );
            extrapolation_bc=zero(eltype(receiver.coil_sens)),
        ) for coil in axes(receiver.coil_sens, 4)
    ]
    return (; interpolated_dims, singleton_dims, interpolations)
end

function get_sens(receiver::ArbitraryCoilSens, x, y, z)
    sens = similar(x, eltype(receiver.coil_sens), length(x), get_n_coils(receiver))
    return get_sens!(sens, receiver, x, y, z)
end

function get_sens!(sens, receiver::ArbitraryCoilSens, x, y, z)
    return get_sens!(sens, receiver, x, y, z, coil_interpolators(receiver))
end

function get_sens!(sens, receiver::ArbitraryCoilSens, x, y, z, interpolators)
    positions = (x, y, z)
    if isempty(interpolators.interpolated_dims)
        sens .= reshape(receiver.coil_sens[1, 1, 1, :], 1, :)
    else
        queries = Tuple(positions[dim] for dim in interpolators.interpolated_dims)
        for coil in axes(receiver.coil_sens, 4)
            sens[:, coil] .= interpolators.interpolations[coil].(queries...)
        end
    end
    for dim in interpolators.singleton_dims
        coordinate = only((receiver.x, receiver.y, receiver.z)[dim])
        for coil in axes(sens, 2), spin in eachindex(positions[dim])
            if positions[dim][spin] != coordinate
                sens[spin, coil] = zero(eltype(sens))
            end
        end
    end
    return sens
end

struct PathInterpolation{P,K}
    coordinates::P
    knots::K
end

struct PathCoordinates{X,D,P,S}
    position::X
    displacement::D
    points::P
    spins::S
end

@adapt_structure PathCoordinates

motion_enabled(::PathInterpolation) = Val(true)

prealloc_motion_coordinates(motion::Motion, backend::KA.GPU, obj, max_block_length) =
    prealloc_motion_coordinates(motion.action, backend, obj, max_block_length)
prealloc_motion_coordinates(::KomaMRIBase.AbstractAction, backend::KA.GPU, obj, max_block_length) =
    dense_motion_coordinates(backend, obj, max_block_length)
function prealloc_motion_coordinates(
    action::Union{Path,FlowPath}, backend::KA.GPU, obj, max_block_length,
)
    1 < size(action.dx, 2) == size(action.dy, 2) == size(action.dz, 2) ||
        return dense_motion_coordinates(backend, obj, max_block_length)
    T = eltype(obj.x)
    knots = gpu(collect(KomaMRIBase.trajectory_knots(size(action.dx, 2), T)), backend)
    points = similar(obj.x, Tuple{Int32,T}, max_block_length)
    spins = KomaMRIBase.expand(obj.motion.spins, length(obj.x)).range
    coordinates = map((position, displacement) -> PathCoordinates(position, displacement, points, spins),
        (obj.x, obj.y, obj.z), (action.dx, action.dy, action.dz))
    return PathInterpolation(coordinates, knots)
end

function spin_coordinates!(interpolation::PathInterpolation, motion, x, y, z, t)
    points = @view interpolation.coordinates[1].points[1:length(t)]
    ts = KomaMRIBase.unit_time(t, motion.time)
    points .= KomaMRIBase.trajectory_point.(Ref(interpolation.knots), vec(ts))
    return interpolation.coordinates
end

@inline function Base.getindex(coordinates::PathCoordinates, i, t)
    @inbounds position = coordinates.position[i]
    spins = coordinates.spins
    i in spins || return position
    spin = size(coordinates.displacement, 1) == 1 ? 1 : (Int(i) - first(spins)) ÷ step(spins) + 1
    @inbounds return position + KomaMRIBase.trajectory_sample(coordinates.displacement, spin, coordinates.points[t])
end

@inline get_spin_coordinates(x::PathCoordinates, y::PathCoordinates, z::PathCoordinates, i, t) =
    (x[i, t], y[i, t], z[i, t])

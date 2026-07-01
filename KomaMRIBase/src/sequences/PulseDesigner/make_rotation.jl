"""
    seq = build_rotation(args...; sys=Scanner())

Return a one-block `Sequence` with a Pulseq-style rotation extension. See
`make_rotation` for rotation arguments.

# Keywords
- `sys=Scanner()`: Scanner defaults and raster times.

# Returns
- `seq`: Sequence containing the rotation block.
"""
function build_rotation(args...; sys=Scanner())
    seq = Sequence(sys)
    addblock!(seq, make_rotation(args...))
    return seq
end

"""
    event = make_rotation(phi)
    event = make_rotation(phi, theta)
    event = make_rotation(axis, angle)
    event = make_rotation(quaternion)
    event = make_rotation(matrix)

Return a Pulseq-style 3D rotation extension.

# Arguments
- `phi`: Rotation about z. [`rad`]
- `theta`: Rotation about y before `phi`. [`rad`]
- `axis`: 3-vector rotation axis.
- `angle`: Rotation about `axis`. [`rad`]
- `quaternion`: 4-vector `[q0, qx, qy, qz]`.
- `matrix`: 3x3 rotation matrix.

# Returns
- `event`: `QuaternionRot` extension event.
"""
make_rotation(phi) = make_rotation(phi, 0.0)

function make_rotation(phi, theta)
    phi = to_SI(phi, SIUnitsDefault())
    theta = to_SI(theta, SIUnitsDefault())
    check_rotation_angle(phi, -π, 2π, "phi"; include_hi=false)
    check_rotation_angle(theta, -π, π, "theta")
    q_y = QuaternionRot(cos(theta / 2), 0, sin(theta / 2), 0)
    q_z = QuaternionRot(cos(phi / 2), 0, 0, sin(phi / 2))
    return q_z * q_y
end

function make_rotation(axis::AbstractVector, angle)
    angle = to_SI(angle, SIUnitsDefault())
    length(axis) == 3 || error("Rotation axis must have length 3.")
    check_rotation_angle(abs(angle), 0, π, "angle")
    norm_axis = sqrt(sum(abs2, axis))
    norm_axis > 0 || error("Rotation axis must be non-zero.")
    axis = axis ./ norm_axis
    return QuaternionRot(cos(angle / 2), (sin(angle / 2) .* axis)...)
end

function make_rotation(quaternion::AbstractVector)
    length(quaternion) == 4 || error("Rotation quaternion must have length 4.")
    qnorm = sqrt(sum(abs2, quaternion))
    qnorm > 0 || error("Rotation quaternion must be non-zero.")
    return QuaternionRot((quaternion ./ qnorm)...)
end

make_rotation(matrix::AbstractMatrix) = QuaternionRot(matrix)
make_rotation(q::QuaternionRot) = q

function check_rotation_angle(angle, lo, hi, name; include_hi=true)
    ok = include_hi ? lo <= angle <= hi : lo <= angle < hi
    interval = include_hi ? "[$lo, $hi]" : "[$lo, $hi)"
    ok || error("Rotation angle `$name` must be in $interval radians.")
    return nothing
end

function Base.:*(q1::QuaternionRot, q2::QuaternionRot)
    q0 = q1.q0 * q2.q0 - q1.qx * q2.qx - q1.qy * q2.qy - q1.qz * q2.qz
    qx = q1.q0 * q2.qx + q1.qx * q2.q0 + q1.qy * q2.qz - q1.qz * q2.qy
    qy = q1.q0 * q2.qy - q1.qx * q2.qz + q1.qy * q2.q0 + q1.qz * q2.qx
    qz = q1.q0 * q2.qz + q1.qx * q2.qy - q1.qy * q2.qx + q1.qz * q2.q0
    return make_rotation([q0, qx, qy, qz])
end

"""
    QuaternionRot(q0, qx, qy, qz)
    QuaternionRot(R)

Pulseq `ROTATIONS` extension stored as a quaternion. `QuaternionRot(R)` converts
a 3x3 rotation matrix using the same convention as Pulseq MATLAB's
`mr.aux.quat.fromRotMat`.
"""
struct QuaternionRot <: Extension
    q0::Float64
    qx::Float64
    qy::Float64
    qz::Float64
end

function QuaternionRot(R::AbstractMatrix{T}) where {T<:Real}
    size(R) == (3, 3) || throw(DimensionMismatch("Rotation matrix must be 3x3."))
    R = copyto!(similar(R, typeof(float(one(T)))), R)
    all(iszero, R) && throw(ArgumentError("Empty matrix provided in place of a rotation matrix."))

    qs = sqrt(max(0, R[1, 1] + R[2, 2] + R[3, 3] + 1)) / 2
    if abs(qs) <= eps(eltype(R))
        q0 = zero(qs)
        sgn_R23 = -R[2, 3] < 0 ? -1 : 1
        sgn_R13 = -R[1, 3] < 0 ? -1 : 1
        sgn_R12 = -R[1, 2] < 0 ? -1 : 1
        qx = sqrt(max(0, -(R[2, 2] + R[3, 3]) / 2)) * sgn_R23
        qy = sqrt(max(0, -(R[1, 1] + R[3, 3]) / 2)) * sgn_R13
        qz = sqrt(max(0, -(R[1, 1] + R[2, 2]) / 2)) * sgn_R12
    else
        q0 = qs
        qx = (R[3, 2] - R[2, 3]) / (4 * qs)
        qy = (R[1, 3] - R[3, 1]) / (4 * qs)
        qz = (R[2, 1] - R[1, 2]) / (4 * qs)
    end

    qnorm = sqrt(abs2(q0) + abs2(qx) + abs2(qy) + abs2(qz))
    qnorm > 0 || throw(ArgumentError("Rotation matrix cannot be converted to a non-zero quaternion."))
    return QuaternionRot(q0 / qnorm, qx / qnorm, qy / qnorm, qz / qnorm)
end

get_scale(::Type{QuaternionRot}) = [1.0 1.0 1.0 1.0]
get_pulseq_format(::Type{QuaternionRot}) = "%f %f %f %f"
get_EXT_type_from_symbol(::Val{:ROTATIONS}) = QuaternionRot
get_symbol_from_EXT_type(::Type{QuaternionRot}) = "ROTATIONS"
extension_type_header(::Type{QuaternionRot}) = "# Extension specification for rotations:\n# id q0 qx qy qz\n"

function rotation_matrix(q::QuaternionRot)
    qnorm = sqrt(abs2(q.q0) + abs2(q.qx) + abs2(q.qy) + abs2(q.qz))
    qnorm > 0 || throw(ArgumentError("Rotation quaternion must have non-zero norm."))
    q0, qx, qy, qz = q.q0 / qnorm, q.qx / qnorm, q.qy / qnorm, q.qz / qnorm
    return [
        1 - 2(qy^2 + qz^2)  2(qx * qy - q0 * qz)  2(qx * qz + q0 * qy)
        2(qx * qy + q0 * qz)  1 - 2(qx^2 + qz^2)  2(qy * qz - q0 * qx)
        2(qx * qz - q0 * qy)  2(qy * qz + q0 * qx)  1 - 2(qx^2 + qy^2)
    ]
end

Base.transpose(q::QuaternionRot) = QuaternionRot(q.q0, -q.qx, -q.qy, -q.qz)
Base.adjoint(q::QuaternionRot) = transpose(q)

Base.:*(q::QuaternionRot, grads::AbstractVector{<:Grad}) = rotation_matrix(q) * grads
